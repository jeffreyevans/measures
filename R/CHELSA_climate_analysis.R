# Converts temperature from Kelvin to Fahrenheit 
# Preciptattion is in kg m^-2 (mm) so for monthly values;
#   kg/m2 = mm; 3600 seconds/hour * 24 hours/day * 30 days/month = 2,592,000
# Precipitation rate could be specified as the mass flow rate of liquid or solid water 
# across a horizontal plane per unit time: Mw in kg m-2 s-1. Water density is a function 
# of temperature but that can be ignored in this context; then the volume flow rate, or 
# precipitation rate, becomes R = Mw/pw in m s-1 or, more conveniently, in units of 
# mm hr-1 or mm day-1. Precipitation rate is the depth to which a flat horizontal surface 
# would have been covered per unit time if no water were lost by run-off, evaporation, 
# or percolation. 
suppressMessages(
  lapply(c("sf", "spatialEco", "terra", "duckdb", "rts", 
           "climetrics", "classInt", "RColorBrewer"), require, 
		   character.only = TRUE))
sf_use_s2(FALSE)

( cty <- c("bhutan", "costa_rica", "canada", "colombia", "peru", "brazil")[5] )

norm.years <- c(1991:2019)
subset.norms = c(TRUE,FALSE)[2]
clr <- colorRampPalette(c("Maroon", "Khaki", "yellow", "lightblue", "blue", "MidnightBlue"))

# Desginated historic/current climate periods 
tp1 = ifelse(subset.norms == TRUE, '1991/2000', '1979/1990')  
tp2 = '2010/2019' 

root = file.path("C:/evans/PFP", cty)
setwd(root)
  dat.dir = file.path(root, "data")
  clim.dir = file.path(root, "climate")
  clim.data.dir = file.path("E:/climate/chelsa") 
dir.create(clim.dir, showWarnings = FALSE)

#********************************************************
# Functions
k.convert <- function(k, m=c("fahrenheit", "celsius")) { 
  switch(m[1], 
    "fahrenheit" = { k = (( k / 10 - 273.15) * 1.8) + 32 }, 
    "celsius" = { k = (k/10) - 273.15 }) 	
  return(k)   
}

climateTrend <- function(r, th = 10) {
  y <- terra::as.data.frame(r, cells = TRUE)
    cell.idx <- y$cell
	  y <- y[,-1]
        didx <- which(rowSums(is.na(y))!= ncol(y))    
        y <- t(y[didx, ])
  N <- apply(y, 2, function(x) sum(!is.na(x)))
    ind <- which(N >= th)
      y <- y[,ind] 
        N <- apply(y, 2, function(x) sum(!is.na(x)))
    x <- matrix(nrow = terra::nlyr(r), ncol = ncol(y))
  x[] <- 1:terra::nlyr(r)
  x1 <- y
    x1[!is.na(x1)] <- 1
      x <- x*x1
  sx <- apply(x, 2, sum, na.rm = T)
    sy <- apply(y, 2, sum, na.rm = T)
      sxx <- apply(x, 2, function(x) sum(x^2, na.rm = T))
      syy <- apply(y, 2, function(x) sum(x^2, na.rm = T))
    xy <- x*y
  sxy <- apply(xy, 2, sum, na.rm = T)
  slope <- (N * sxy - (sx * sy)) / (N * sxx-sx^2)
    sres <- (N * syy - sy^2 - slope^2 * (N * sxx - sx^2)) / (N * (N - 2))
      SE <- suppressWarnings(sqrt((N * sres) / (N * sxx - sx^2)))
    Test <- slope / SE
  p <- mapply(function(x,y) (2 * pt(abs(x), df = y - 2, lower.tail = FALSE)), x = Test, y = N)
    terra::time(r) <- NULL
      slpTrends <- sigTrends <- seTrends <- terra::rast(r[[1]])
	    slpTrends[] <- NA
	      sigTrends[] <- NA
	        seTrends[] <- NA
            slpTrends[cell.idx[ind]] <- slope
          seTrends[cell.idx[ind]] <- SE
        sigTrends[cell.idx[ind]] <- p
      output <- c(slpTrends, seTrends, sigTrends)
    names(output) <- c("slpTrends", "seTrends", "sigTrends")
  return(output)
}

timeShift <- function(r, yr0, yr1, yr2, m, trend.method = c("schoeman", "taheri")){
  m1 <- ((yr1 - yr0) * 12) + m
    m2 <- ((yr2 - yr0) * 12) + m
      r1 <- r[[seq(m1, m2, by = 12)]]
    switch(trend.method[1], 
      "schoeman" = { trend <- climateTrend(r1, th = 10)[[1]] },  
      "taheri" = { trend <- climetrics::temporalTrend(rts::rts(r1, 
                              time = as.Date(names(r1))))[[1]] } )
  b <- ifelse((m - 1) == 0, 12, (m - 1))
    m1 <- ((yr1-yr0)*12)+b
    m2 <- ((yr2-yr0)*12)+b
  x2 <- r[[seq(m1, m2, by = 12)]]
  b <- ifelse((m + 1) == 13, 1, (m + 1))
    m1 <- ((yr1 - yr0) * 12) + b
    m2 <- ((yr2 - yr0) * 12) + b
  x3 <- r[[seq(m1, m2, by = 12)]]
  x4 <- terra::mean((x3 - x2)/2, na.rm = TRUE)
	if(trend.method[1] == "taheri") { 
	  sShift <- ( (trend*0.1) / x4) * (3652.5 / 12)
    } else { 
	  sShift <- (trend / x4) * (3652.5 / 12)
	}
	sShift[sShift == Inf | sShift == -Inf] <- NA
        timing <- c(trend, x4, sShift)
          names(timing) <- c("trend", "rate", "shift")
  return(timing)
}

#********************************************************
# If needed, create climate metric
# for temp Kelvin is converted to Fahrenheit 
bdy <- st_read(file.path(dat.dir, paste0(cty, ".gpkg")), "boundary")    
  e <- ext(st_transform(bdy, 4326)) 
    ref <- rast(ext(bdy), resolution = 1000, crs = crs(vect(bdy)))
      ref[] <- 1
	    ref <- mask(ref, bdy)

# Min temp
rname <- file.path(clim.dir, paste0("temp_min.tif"))
  if(!file.exists(rname)) {
    f <- list.files(file.path(clim.data.dir, "tasmin"), "tif$", full.names=TRUE)
    dates <- as.Date(unlist(lapply(strsplit(basename(f), "_"), \(i) {
               paste(i[4], i[3], "01", sep="/")})))
      date.idx <- sort.int(dates, index.return = TRUE)$ix
        dates <- dates[date.idx]
    clim <- rast(f[date.idx])
      terra::time(clim, tstep="months") <- dates
        names(clim) <- dates  
          clim <- crop(clim, e)
            clim <- mask(project(clim, ref, method="bilinear"), ref)
              clim <- app(clim, k.convert)
    writeRaster(clim, rname)
  }
  
# Max temp
rname <- file.path(clim.dir, paste0("temp_max.tif"))
  if(!file.exists(rname)) {
    f <- list.files(file.path(clim.data.dir, "tasmax"), "tif$", full.names=TRUE)
    dates <- as.Date(unlist(lapply(strsplit(basename(f), "_"), \(i) {
               paste(i[4], i[3], "01", sep="/")})))
      date.idx <- sort.int(dates, index.return = TRUE)$ix
        dates <- dates[date.idx]
    clim <- rast(f[date.idx])
      terra::time(clim, tstep="months") <- dates
        names(clim) <- dates  
          clim <- crop(clim, e)
            clim <- mask(project(clim, ref, method="bilinear"), ref)
              clim <- app(clim, k.convert)
    writeRaster(clim, rname)
  } 

# Mean temp
rname <- file.path(clim.dir, paste0("temp_mean.tif"))
  if(!file.exists(rname)) {
    f <- list.files(file.path(clim.data.dir, "tas"), "tif$", full.names=TRUE)
    dates <- as.Date(unlist(lapply(strsplit(basename(f), "_"), \(i) {
               paste(i[4], i[3], "01", sep="/")})))
      date.idx <- sort.int(dates, index.return = TRUE)$ix
        dates <- dates[date.idx]
    clim <- rast(f[date.idx])
      terra::time(clim, tstep="months") <- dates
        names(clim) <- dates  
          clim <- crop(clim, e)
            clim <- mask(project(clim, ref, method="bilinear"), ref)
              clim <- app(clim, k.convert)
    writeRaster(clim, rname)
  } 

# Precip
rname <- file.path(clim.dir, paste0("precip.tif"))
  if(!file.exists(rname)) {
    f <- list.files(file.path(clim.data.dir, "pr"), "tif$", full.names=TRUE)
    dates <- as.Date(unlist(lapply(strsplit(basename(f), "_"), \(i) {
               paste(i[4], i[3], "01", sep="/")})))
      date.idx <- sort.int(dates, index.return = TRUE)$ix
        dates <- dates[date.idx]
    clim <- rast(f[date.idx])
      terra::time(clim, tstep="months") <- dates
        names(clim) <- dates  
          clim <- crop(clim, e)
            clim <- mask(project(clim, ref, method="bilinear"), ref)
    writeRaster(clim, rname)
  } 

# climate moisture index
rname <- file.path(clim.dir, paste0("cmi.tif"))
  if(!file.exists(rname)) {
    f <- list.files(file.path(clim.data.dir, "cmi"), "tif$", full.names=TRUE)
    dates <- as.Date(unlist(lapply(strsplit(basename(f), "_"), \(i) {
               paste(i[4], i[3], "01", sep="/")})))
      date.idx <- sort.int(dates, index.return = TRUE)$ix
        dates <- dates[date.idx]
    clim <- rast(f[date.idx])
      terra::time(clim, tstep="months") <- dates
        names(clim) <- dates  
          clim <- crop(clim, e)
            clim <- mask(project(clim, ref, method="bilinear"), ref)
    writeRaster(clim, rname)
  } 

# Bioclime metrics (1981-2010 normals)
rname <- file.path(clim.dir, paste0("bioclime.tif"))
  if(!file.exists(rname)) {
    f <- list.files(file.path(clim.data.dir, "bioclime"), "tif$", full.names=TRUE)
      vnames <- unlist(lapply(strsplit(basename(f), "_"), \(i) { 
        if(length(i) == 4) {
          return(i[2]) 
        } else if(length(i) == 5) {
          return(paste0(i[2],"_", i[3]))
        } else if(length(i) == 6) {
          return(paste0(i[2],"_", i[3], "_", i[4]))
        }		
      }))
	  clim <- rast(f)
        names(clim) <- vnames
          clim <- crop(clim, e)
            clim <- mask(project(clim, ref, method="bilinear"), ref)
    writeRaster(clim, rname)
  } 

#********************************************************
# read rasters and subset to climate normal period 1991-2019
tmin <- rast(file.path(clim.dir, "temp_min.tif"))
  if(subset.norms) tmin <- tmin[[grep(paste(norm.years, collapse="|"),names(tmin))]]
    dates <- as.Date(names(tmin))
      terra::time(tmin, tstep="months") <- as.Date(names(tmin))	

tmax <- rast(file.path(clim.dir, "temp_max.tif"))
  if(subset.norms) tmax <- tmax[[grep(paste(norm.years, collapse="|"),names(tmax))]]
    dates <- as.Date(names(tmax))
      terra::time(tmax, tstep="months") <- as.Date(names(tmax))	

tmean <- rast(file.path(clim.dir, "temp_mean.tif"))
  if(subset.norms) tmean <- tmean[[grep(paste(norm.years, collapse="|"),names(tmean))]]
    dates <- as.Date(names(tmean))
      terra::time(tmean, tstep="months") <- as.Date(names(tmean))	

precip <- rast(file.path(clim.dir, "precip.tif")) * 0.001
  if(subset.norms) precip <- precip[[grep(paste(norm.years, collapse="|"),names(precip))]]
    dates <- as.Date(names(precip))
      terra::time(precip, tstep="months") <- as.Date(names(precip))	

cmi <- rast(file.path(clim.dir, "cmi.tif")) * 0.001
  if(subset.norms) cmi <- cmi[[grep(paste(norm.years, collapse="|"),names(cmi))]]
    dates <- as.Date(names(cmi))
      terra::time(cmi, tstep="months") <- as.Date(names(cmi))	

# bio <- rast(file.path(clim.dir, "bioclime.tif"))

#********************************************************
# Climate change metrics
# - Standardized Local Anomalies (sed)
# - localExtreme (lce): Changes in probability of local extremes
# - aaClimate (aac): Changes in areas of analogous climates
# - novelClimate (nc): Novel climates
# - daClimate (dac): Changes in distance to analogous climates
# - velocity (ve): Climate change velocity 
# - dVelocity (dVe): Distance-based climate change velocity 
# - gVelocity (gVe): Gradiant-based climate change velocity 

# coerce data to rts objects and calculate climate metrics
tmin <- rts(tmin, time = as.Date(names(tmin))) 
tmax <- rts(tmax, time = as.Date(names(tmax)))  
tmean <- rts(tmean, time = as.Date(names(tmean))) 
precip <- rts(precip, time = as.Date(names(precip))) 
cmi <- rts(cmi, time = as.Date(names(cmi)))   

# Shift in timing
if(!file.exists(file.path(clim.dir, "climate_timing.tif"))) {
  timing.tmax <- timeShift(tmax@raster, yr0=1979, yr1=1985, yr2=2019, m=5)
    names(timing.tmax) <- c("tmax_trend", "tmax_rate", "tmax_shift") 
  timing.precip <- timeShift(precip@raster, yr0=1979, yr1=1985, yr2=2019, m=5)
    names(timing.precip) <- c("precip_trend", "precip_rate", "precip_shift") 
  writeRaster(timing, file.path(clim.dir, "climate_timing.tif"), overwrite=TRUE)
} else {
  timing <- rast(file.path(clim.dir, "climate_timing.tif"))
}

# Climate trends (min, mean, max temp, precip and cmi)
if(!file.exists(file.path(clim.dir, "climate_trends.tif"))) {
  trend.precip <- timing.precip[[1]] 
    names(trend.precip) <- c("precip_trend") 
  trend.tmax <- timing.tmax[[1]] 
    names(trend.tmax) <- c("tmax_trend") 
  trend.tmin <- climateTrend(tmin@raster)[[1]] * 10
    names(trend.tmin) <- c("tmin_trend") 
  trend.tmean <- climateTrend(tmean@raster)[[1]] * 10
    names(trend.tmean) <- c("tmean_trend")   
  trend.cmi <- climateTrend(cmi@raster)[[1]] * 10
    names(trend.cmi) <- c("cmi_trend") 
  trends <- c(trend.precip, trend.tmin, trend.tmean, trend.tmax, trend.cmi)
  timing <- c(timing.tmax[[-1]], timing.precip[[-1]])
  writeRaster(trends, file.path(clim.dir, "climate_trends.tif"), overwrite=TRUE)
} else {
  trends <- rast(file.path(clim.dir, "climate_trends.tif"))
}

if(!file.exists(file.path(clim.dir, "climate_metrics.tif"))) {
  lce <- localExtreme(precip, tmax, t1 = tp1, t2 = tp2, extreme = c(0.95,0.05))  
  se <- sed(precip, tmin, tmax, t1=tp1, t2=tp2) 
  ve <- velocity(x1=precip, x2=tmax, t1=tp1, t2=tp2)
  dv <- dVelocity(precip, tmin, tmax, t1=tp1, t2=tp2) 	   
  clim.metrics <- c(lce, se, ve, dv)
    names(clim.metrics) <- c("localExtreme", "anomalies", "velocity", "dVelocity")
  writeRaster(clim.metrics, file.path(clim.dir, "climate_metrics.tif"), overwrite=TRUE)
} else {
  clim.metrics <- rast(file.path(clim.dir, "climate_metrics.tif")) 
}

pdf(file.path(clim.dir, "climate_plots.pdf"), height=10, width=10)
  cl <- colorRampPalette(c("MidnightBlue","Turquoise","lightblue","gray","khaki","orange","red"))
  # probability of extreme climate events
  if("localExtreme" %in% names(clim.metrics))
    plot(clim.metrics[["localExtreme"]], main='probability of extreme climate events', 
         col=rev(clr(100)), maxcell=5000000, smooth=TRUE)
   
  # Standardized Local Anomalies
  if("anomalies" %in% names(clim.metrics))  
    plot(clim.metrics[["anomalies"]], main='Standardized Local Anomalies', 
         col=rev(clr(100)), maxcell=5000000, smooth=TRUE)
  
  # Climate change Velocity
  if("velocity" %in% names(clim.metrics))  
    plot(clim.metrics[["velocity"]], col=cl(100), main='Velocity of Climate Change', 
         maxcell=5000000, smooth=TRUE)
  
  # Distance-based Velocity
  if("dVelocity" %in% names(clim.metrics))   
    plot(clim.metrics[["dVelocity"]], col=cl(100), main='Distance-based Velocity of Climate Change',
         maxcell=5000000, smooth=TRUE)
  
  # precip trend
  if("precip_trend" %in% names(trends))  
    plot(trends[["precip_trend"]], col=rev(cl(100)), main='precipitation trend 1979-2018', 
         maxcell=5000000, smooth=TRUE)

  # moisture index trend
  if("cmi_trend" %in% names(trends))
    plot(trends[["cmi_trend"]], col=rev(cl(100)), main='moisture index trend 1979-2018', 
         maxcell=5000000, smooth=TRUE)	 
		 
  # min temp trend
  if("tmin_trend" %in% names(trends))    
    plot(trends[["tmin_trend"]], col=cl(100), main='minimum temperature trend 1979-2018', 
         maxcell=5000000, smooth=TRUE)

  # mean temp trend
  if("tmean_trend" %in% names(trends))      
    plot(trends[["tmean_trend"]], col=cl(100), main='mean temperature trend 1979-2018', 
         maxcell=5000000, smooth=TRUE)

  # max temp trend
  if("tmax_trend" %in% names(trends))      
    plot(trends[["tmax_trend"]], col=cl(100), main='maximum temperature trend 1979-2018', 
         maxcell=5000000, smooth=TRUE)

  # precip change in timing rate
  if("precip_rate" %in% names(timing))      
    plot(timing[["precip_rate"]], col=cl(100), main='Rate of change in timing for precipitation', 
         maxcell=5000000, smooth=TRUE)

  # max temp change in timing rate
  if("tmax_rate" %in% names(timing))      
    plot(timing[["tmax_rate"]], col=cl(100), main='Rate of change in timing for maximum temperature', 
         maxcell=5000000, smooth=TRUE)
	   
dev.off()

tmpFiles(current=TRUE, orphan=TRUE, old=TRUE, remove=TRUE)
