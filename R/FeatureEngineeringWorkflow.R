# Feature engineering to produce design matrix covariates
#   that represent characteristics of the timeseries
# https://cran.r-project.org/web/packages/tsfeatures/vignettes/tsfeatures.html
suppressMessages(
  lapply(c("sf", "spatialEco", "terra", "dplyr", "lubridate", 
         "grf", "rts", "exactextractr", "tsfeatures"), 
		 require, character.only = TRUE))

#*************************************************
# required functions
# x <- as.numeric(intervention[1,][8:ncol(intervention)])
heterogeneity <- function(x) {
    x.whitened <- na.contiguous(ar(x)$resid)
    x.archtest <- arch_stat(as.matrix(x.whitened))
    LBstat <- sum(acf(x.whitened^2, lag.max = 12L, plot = FALSE)$acf[-1L]^2)
    garch.fit <- suppressWarnings(tseries::garch(x.whitened, trace = FALSE))
    garch.fit.std <- residuals(garch.fit)
    x.garch.archtest <- arch_stat(garch.fit.std)
    LBstat2 <- NA
    try(LBstat2 <- sum(acf(na.contiguous(garch.fit.std^2), lag.max = 12L, 
        plot = FALSE)$acf[-1L]^2), silent = TRUE)
    output <- c(arch_acf = LBstat, garch_acf = LBstat2, arch_r2 = unname(x.archtest), 
        garch_r2 = unname(x.garch.archtest))
    return(output)
}

process.ts <- function(j, start = c(2014, 1), f=36) {
  ts.features <- function(x) {
    tsf <- c(entropy(x), stability(x), lumpiness(x),
      hurst(x), max_level_shift(x), max_var_shift(x),
      max_kl_shift(x), crossing_points(x),
      flat_spots(x), stl_features(x)[c(3,4,6,7,8)],
      acf_features(x), pacf_features(x), arch_stat(x),
      heterogeneity(x), firstmin_ac(x))
    return(tsf)
  } 
  j <- ts(j, start=start, frequency=f)
      ts.metrics <- ts.features(j)  
    names(ts.metrics)[34] <- "facf"  
  return(ts.metrics)
} 

#*************************************************
# Feature Engineering function
feature.engineering <- function(y) {


#*************************************************
# Set environment and model parameters 
metrics <- c("lai", "fcov")[1]

# dry 12 and 1:5 and wet 6:11
season = c("none", "dry", "wet")[2] 

root = "C:/evans/Colombia/waterfunds"
dat.dir = file.path(root, "data")
mdl.dir <- file.path(root, "model", "causal")
setwd(root)

#**************************************
# read point data with control and treatment
#   intervention n=831 (%forest > 20% n=293) 
#   control n=1641 (%forest > 20% n=915) 
#   protected areas n=7076 (%forest > 20% n=2462)  
# If lai, subset forest > 20%   
cells <- st_read(file.path(dat.dir, "model_300m_data.gpkg"))
  st_geometry(cells) <- "geometry"
  dat <- na.omit(cells)  
    dat$intervention <- ifelse(dat$protected.area != "no", 
	                      "Protected Area", dat$intervention)  
      if(metrics == "lai") {
        dat <- dat[dat$pct.forest >= 0.20,] 
      }

#*************************************************
# read raster time-series and parse dates
# extract time series to model data, fill NA's in 
# the timeseries and then remove > 6 missing values
metric <- rast(file.path(dat.dir, paste0(toupper(metrics), "_300m_trend", ".tif")))
  d <- unlist(lapply(strsplit(names(metric), "_"), function(x) x[2]))
  dates <- as.Date(unlist(lapply(d, function(j) {
      paste( substr(j,1,4), substr(j,5,6), substr(j,7,8), sep="-" )
    })))
	if(metrics == "lai") {
	  if(min(global(metric, fun="min", na.rm=TRUE)[,1]) < 0)
        metric[metric < 0] <- 0
	  if(max(global(metric, fun="max", na.rm=TRUE)[,1]) > 9)
	    metric[metric > 9] <- 9
    } else if(metrics == "fcov") {
	  if(min(global(metric, fun="min", na.rm=TRUE)[,1]) < 0)
        metric[metric < 0] <- 0
	  if(max(global(metric, fun="max", na.rm=TRUE)[,1]) > 1)
	    metric[metric > 1] <- 1
    }	

#*************************************************
# create timeseries data.frame object
tdat <- terra::extract(metric, vect(dat))[,-1]  
  names(tdat) <- gsub("LAI", toupper(metrics), names(tdat))
    na.idx <- as.data.frame(which(is.na(tdat), arr.ind = TRUE))
	if(nrow(na.idx) > 0) {
      na.cts <- as.data.frame(table(na.idx$row))
	      na.cts <- na.cts[na.cts$Freq > 4,]
          na.idx <- as.numeric(as.character(na.cts[which(na.cts$Freq > 10),]$Var1)) 
        if(length(na.idx) > 0) {
	      dat <- dat[-na.idx,]
	      tdat <- tdat[-na.idx,]
	    }
	}

#*************************************************
# Time-series Feature Engineering on controls
if(!file.exists(file.path(mdl.dir, paste0("Control_TSFE_", metrics, ".rds")))){
  cf.idx <- which(dat$intervention == "no")
  ctse <- setNames(vector("list", length(cf.idx)), dat$cell[cf.idx])
    for(i in 1:length(cf.idx)) {
	  v <- cf.idx[i]
      tryCatch({
        fe <- process.ts(na.omit(as.numeric(tdat[v,])))
          }, error=function(e){ fe <- rep(NA,34) })
      ctse[[v]] <- fe
    }	
  ctse <- do.call(rbind, ctse)
    ctse <- data.frame(cell=dat$cell[cf.idx], ctse)  
      saveRDS(ctse, file = file.path(mdl.dir, 
	          paste0("Control_TSFE_", metrics, ".rds")))
} 

#*************************************************
# Time-series Feature Engineering on interventions    
if(!file.exists(file.path(mdl.dir, paste0(j, "_Intervention_TSFE_", metrics, ".rds")))){
itse <- setNames(vector("list", length(int.idx)), dat$cell[int.idx])
  for(i in 1:length(int.idx)) {
  v <- int.idx[i]
    tryCatch({
      fe <- process.ts(na.omit(as.numeric(tdat[v,])))
        }, error=function(e){ fe <- rep(NA,34) })
    itse[[v]] <- fe	
  }	 	
  itse <- do.call(rbind, itse)
    itse <- data.frame(cell=dat$cell[int.idx], itse)
      saveRDS(itse, file = file.path(mdl.dir, paste0("Intervention_TSFE_", 
	        metrics, ".rds")))
}
