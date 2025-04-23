# process time series filling NA's, optionally smoothing 
# and decomposing periodicity to derive a trend. This uses
# a linear trend method or multi-model Bayesian averaging 
# via the BEAST method
#
# raster_300m_filled  - Imputed NA's and smoothed. The default
#   is a Kalman filter. 
# raster_300m_trend - Decomposed trend. In the case of
#   beast (Bayesian), the trend in non-linear 
#
library(terra)
library(Rbeast)
library(spatialEco)
library(imputeTS)
library(bfast)

idx = 2  # indicates country
  country <- c("bhutan", "costa_rica", "canada", "colombia", "peru", "brazil")[idx]

metrics <- c("lai", "fcov")[2]
mres = c("250m", "300m", "500m")[3]
  if(mres == "250m") metric.type = "lai"
  if(mres == "500m") metric.type = "fcov"

out.file <- paste0(metrics, "_", mres, "_trend.tif") 
in.file <- paste0(metrics, "_", mres, "_raw.tif")  

root = file.path("C:/evans/PFP/timeseries", country)
dat.dir = file.path(root, paste0(metrics, mres))
setwd(root)

model = c("moving.averages", "Bayesian", "bfast")[1]
smoothing = c("kalman", "lowess")[2]
smooth.data = FALSE
mobs = 36                 # minimum number of observsations in timeseries
gap.size = c(Inf, 20)[1]  # maxmum gap size of missing data, defaults to none   

# apply forest mask, p is fraction that timeseires is under threshold
forest.mask = c(TRUE, FALSE)[1]
forest.threshold = c(2, 2.5, 1.5, 2.5, 2.5)[idx]
p = 0.30

#*********************************************************
# Functions
detrend <- function(x, fdate = c(2014, 1), freq = 36, min.obs = 10, 
                    rm.na = FALSE, start.date = "2014-01-10", 
                    periodicity = c("additive", "multiplicative", "beast", "bfast")) {
  if(!any(periodicity[1] == c("additive", "multiplicative", "beast", "bfast")))
    stop("Not a valid periodicity method")
  if(length(x[!is.na(x)]) < min.obs) {
    x.trend <- rep(NA, length(x))
  } else {
	if(any(periodicity[1] == c("additive", "multiplicative"))) {
	  x.ts = stats::ts(na.omit(x), start = fdate, frequency = freq) 
	    x.trend <- as.numeric(stats::decompose(x.ts, type = periodicity[1])$trend)
	} else if(periodicity[1] == "befast") {
	  x.ts = ts(na.omit(x), start = fdate, frequency = freq) 
        d <- bfast::bfastpp(x.ts, stl = "both", lag = 1:2)
          x.trend <- stats::lm(response ~ lag, data = d)$fitted.values
	} else if(periodicity[1] == "beast") {
	  x.trend <- Rbeast::beast(as.numeric(x), period=9, deltat=3/12, 
	                           start=as.Date(start.date), 
							   mcmc.samples = 1000,
							   print.options  = FALSE,
                               print.progress = FALSE, 
							   quiet = TRUE)
		x.trend <- x.trend$trend$Y
	}
	if(!rm.na) {
	  if(length(x) > length(x.trend))  
	    x.trend <- spatialEco::insert.values(x.trend, NA, which(is.na(x))) 
    }	  
  }  
  return( x.trend )
}

fill.na <- function(x, min.obs = mobs) {
  if(length(x[!is.na(x)]) >= min.obs) {
  x <- na_kalman(x, model = "StructTS", smooth = smooth.data, 
                 nit = -1, maxgap = gap.size)
  } else {
    message("too few observsations, returning NA's")
	x <- rep(NA, length(x))
  }
  return(x)  
}

#***************************************
# read raster time-series and parse dates
if(file.exists(in.file)) { 
  metric <- rast(in.file)
    if(mres == "300m") {
    d <- unlist(lapply(strsplit(names(metric), "_"), function(x) x[2]))
      dates <- as.Date(unlist(lapply(d, function(j) {
        paste( substr(j,1,4), substr(j,5,6), substr(j,7,8), sep="-" )
      })))
	} else {
      dates <- as.Date(names(metric))
    }	  
} else {
  warning("Raw timeseries ", in.file, " does not exist")
}


## calculate percentate of timeseries data
#  pdat <- app(metric, fun=\(i) length(i[!is.na(i)]) / length(i) )
#    writeRaster(pdat, paste0(toupper(metrics), "_300m_frac_data", ".tif"),
#              overwrite=TRUE)
# #pidx <- ifel(pdat > 0.20, 1, NA)
# #pidx <- which(fdat > 0.10)

if(pct.mask) {
  cat("**** Masking to timeseries with > 30% data **** ", "\n")
  pct <- rast(paste0(toupper(metrics), "_300m_", "frac_data", ".tif"))
    pct <- ifel(pct >= 0.3, 1, NA)
      metric <- mask(metric, pct)
}

#*********************************************************
#*********************************************************
# Fill NA values in the time-series and write "metric_300m_filled.tif"
# Screens for erronious min-max values and clamps them; fcov 0-1 
# or lai 0-9.5 
if(!file.exists(paste0(toupper(metrics), "_300m_raw", ".tif"))) {
  mm <- terra::global(metric, c("min", "max"), na.rm=TRUE)
    mm <- c(min(mm[,1]), max(mm[,2]))
    if(mm[1] < 0){
      metric[metric < 0] <- 0
    }
    if(metrics == "fcov" & mm[2] > 1) { 
      metric[metric > 1] <- 1
    }
    if(metrics == "lai" & mm[2] > 9.5) {
      metric[metric > 9.5] <- 9
    }	  
  if(smoothing == "lowess") {
    cat("Interpolating NA's", toupper(metrics), "using Lowess regression", "\n")
    metric.fill <- smooth.time.series(metric, f = 0.4, smooth.data = TRUE)
      writeRaster(metric.fill, paste0(toupper(metrics), "_300m_smoothed", ".tif"),
                  overwrite=TRUE)	
  } else if(smoothing == "kalman") {
    cat("Interpolating NA's and smoothing", toupper(metrics), "using a Kalman filter", "\n")  
    #metric <- terra::app(metric, fill.na, cores = 8) 
	tsm <- as.data.frame(metric, cells=TRUE)
      cell.ids <- tsm[,1]
        tsm <- tsm[,-1]
    xts <- apply(tsm, MARGIN=1, FUN=fill.na)
      xts <- as.data.frame(t(xts))
	mm <- c(min(xts,na.rm=TRUE), max(xts,na.rm=TRUE)) 
	  if(mm[1] < 0){ xts[xts < 0] <- 0 } 
      if(metrics == "fcov" & mm[2] > 1) {
		xts[xts > 1] <- 1
      }	  
      if(metrics == "lai" & mm[2] > 9.5) {		
		xts[xts > 9.5] <- 9.5 
	  }	
	metric.fill <- metric[[1]]
      metric.fill[] <- rep(NA, ncell(metric.fill))
        metric.fill <- rep(metric.fill, nlyr(metric))
          names(metric.fill) <- names(metric)  
	metric.fill[cell.ids] <- xts
  writeRaster(metric.fill, paste0(toupper(metrics), "_300m_kmsmoothed", ".tif"),
              overwrite=TRUE)	
  }
} else {
  cat("\n", toupper(metrics), "already filled and smoothed", "\n")
}

#*********************************************************
#*********************************************************
# Decompose periodicity and write trend
# write "metric_300m_trend.tif" 

#*********************************************************
# read required data
if(!file.exists(paste0(toupper(metrics), "_300m_raw", ".tif"))) {
  message(paste0(toupper(metrics), "_300m_raw", ".tif"), " does not exist")  
} else {
  message("reading data ", paste0(toupper(metrics), "_300m_smoothed", ".tif"))
  metric <- rast(paste0(toupper(metrics), "_300m_raw", ".tif"))
    mm <- terra::global(metric, c("min", "max"), na.rm=TRUE)
	  mm <- c(min(mm[,1]), max(mm[,2]))
    if(mm[1] < 0){
      metric[metric < 0] <- 0 
    }
    if(metrics == "fcov" & mm[2] > 1) { 
      metric[metric > 1] <- 1
    }
    if(metrics == "lai" & mm[2] > 9.5) {
      metric[metric > 9.5] <- 9.5
    }	  
#*********************************************************
# using moving averages method
if(model == "moving.averages") {
  if(!file.exists(paste0(toupper(metrics), "_300m_deseason", ".tif"))) {
    cat("\n", "Decomposing periodicity and deriving trend for", toupper(metrics), "\n")
	  cat("Using Moving Averages method", "\n")  
      r.trend <- terra::app(metric, detrend)  
	  na.idx <- global(r.trend, fun="isNA")[,1]
        na.idx <- which(na.idx == ncell(r.trend))
	r.trend <- r.trend[[-na.idx]]
      names(r.trend) <- paste(toupper(metrics), gsub("-", "", dates[-na.idx]), sep="_")	
    terra::writeRaster(r.trend, paste0(toupper(metrics), "_300m_deseason", ".tif"),
                       overwrite=TRUE)
  } else {
    cat("\n", "Periodicity and trend for", toupper(metrics), "already exists", "\n")
  }

#*********************************************************
# using Bayesian (BEAST) method
} else if(model == "Bayesian") {
  if(!file.exists(paste0(toupper(metrics), "_300m_beast_deseason", ".tif"))) {
    cat("\n", "Decomposing periodicity and deriving trend for", toupper(metrics), "\n")
	  cat("Using Bayesian (BEAST) method", "\n")
	  r.trend <- metric
      tsm <- as.data.frame(metric, cells=TRUE)
        rownames(tsm) <- tsm[,1]
          tsm <- tsm[,-1] 
	cat("Running model for n =", nrow(tsm), "\n")
      parms = list()        
        parms$whichDimIsTime = 2
		parms$maxMissingRateAllowed = 0.50 
		parms$period = 1.0
		parms$start=dates[1]
		#parms$time = dates
		if(res == "250m") {
          parms$startTime = dates[1] # c(2020,2,26)
          parms$deltaTime = 1/46 
		} else if(res == "300m") {
		  parms$startTime = dates[1] #c(2014,1,10)
          parms$deltaTime = 1/36         
		}
      o <- beast123(as.matrix(tsm), parms)  
	    mm <- c(min(o$trend$Y[!is.nan(o$trend$Y)], na.rm=TRUE), 
		        max(o$trend$Y[!is.nan(o$trend$Y)], rm=TRUE)) 
	    if(mm[1] < 0){ o$trend$Y[o$trend$Y < 0] <- 0 } 
          if(metrics == "fcov" & mm[2] > 1) {
		    o$trend$Y[o$trend$Y > 1] <- 1
          }	  
          if(metrics == "lai" & mm[2] > 9.5) {		
		    o$trend$Y[o$trend$Y > 9.5] <- 9.0 
	      }	
    r.trend[as.numeric(rownames(tsm))] <- as.data.frame(o$trend$Y)
	  terra::writeRaster(r.trend, paste0(toupper(metrics), "_", res, "_trend", ".tif"),
                         overwrite=TRUE)
  } else {
    cat("\n", "Periodicity and trend for", toupper(metrics), "already exists", "\n")
  }
}

} # end
