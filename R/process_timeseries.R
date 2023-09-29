# process time series filling NA's, optionally smoothing 
# and decomposing periodicity to derive a trend. This uses
# a linear trend method or multi-model Bayesian averaging 
# via the BEAST method
library(terra)
library(Rbeast)
library(spatialEco)
library(imputeTS)
library(bfast)

# setwd("C:/evans/timeseries/CostaRica")
setwd("C:/evans/timeseries/Colombia/CuencaVerde")
# setwd("C:/evans/timeseries/Bhutan")

metrics <- c("lai", "fcov")[2]
model = c("moving.averages", "Bayesian", "bfast")[2]
smoothing = c("kalman", "lowess")[1]

#***************************************
# read raster time-series and parse dates
metric <- rast(paste0(toupper(metrics), "_300m", ".tif"))
  d <- unlist(lapply(strsplit(names(metric), "_"), function(x) x[2]))
  dates <- as.Date(unlist(lapply(d, function(j) {
      paste( substr(j,1,4), substr(j,5,6), substr(j,7,8), sep="-" )
    })))

#*********************************************************
#*********************************************************
# Fill NA values in the time-series and
# write "metric_300m_filled.tif" 
fill.na <- function(x, min.obs = 6) {
  if(length(x[!is.na(x)]) >= min.obs) {
    x <- imputeTS::na_kalman(x)
  }
  return(x)
}

if(!file.exists(paste0(toupper(metrics), "_300m_filled", ".tif"))) {
  if(smoothing == "lowess") {
    cat("filling and smoothing", toupper(metrics), "using Lowess", "\n")
    metric <- smooth.time.series(metric, f = 0.4, smooth.data = FALSE)
  } else if(smoothing == "kalman") {
    cat("filling and smoothing", toupper(metrics), "using Kalman", "\n")  
    metric <- terra::app(metric, fill.na) 
  }
    metric <- ifel(metric < 0, 0, r) 
      if(metrics == "fcov") metric <- ifel(metric > 1, 1, r) 	  
  writeRaster(metric, file.path(root, paste0(toupper(metrics), "_300m_filled", ".tif")),
              overwrite=TRUE)
} else {
  cat("\n", toupper(metrics), "already filled and smoothed", "\n")
}

#*********************************************************
#*********************************************************
# Decompose periodicity and write trend
# write "metric_300m_trend.tif" 
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

#*********************************************************
# using moving averages method
if(model == "moving.averages") {
  if(!file.exists(paste0(toupper(metrics), "_300m_trend", ".tif"))) {
    cat("\n", "Decomposing periodicity and deriving trend for", toupper(metrics), "\n")
	  cat("Using Moving Averages method", "\n")
    metric <- rast(paste0(toupper(metrics), "_300m_filled", ".tif"))
      r.trend <- terra::app(metric, detrend)
        terra::writeRaster(r.trend, paste0(toupper(metrics), "_300m_trend", ".tif"),
                           overwrite=TRUE)
  } else {
    cat("\n", "Periodicity and trend for", toupper(metrics), "already exists", "\n")
  }

#*********************************************************
# using Bayesian (BEAST) method
} else if(model == "Bayesian") {
  if(!file.exists(paste0(toupper(metrics), "_300m_trend", ".tif"))) {
    cat("\n", "Decomposing periodicity and deriving trend for", toupper(metrics), "\n")
	  cat("Using Bayesian (BEAST) method", "\n")
    metric <- rast(paste0(toupper(metrics), "_300m_filled", ".tif"))
      tsm <- as.data.frame(metric, cells=TRUE)
        rownames(tsm) <- tsm[,1]
          tsm <- tsm[,-1] 
	cat("Running model for n =", nrow(tsm), "\n")
      parms = list()        
       parms$whichDimIsTime = 2
       parms$period=9
       parms$deltat=3/12
       parms$start=dates[1] 
      o <- beast123(as.matrix(tsm), parms) 
        metric[as.numeric(rownames(tsm))] <- as.data.frame(o$trend$Y)
      terra::writeRaster(metric, paste0(toupper(metrics), "_300m_trend", ".tif"),
                         overwrite=TRUE)
  } else {
    cat("\n", "Periodicity and trend for", toupper(metrics), "already exists", "\n")
  }
}

#*****************************************************
# pull single pixel time series for testing
## click(metric[[1]], n=1, xy=TRUE, cell=TRUE)
## x=490640.1 y=1124818 cell=475838 
# ( x <- as.numeric(metric[cellFromXY(metric, cbind(490640.1, 1124818))]) )
## x.ts <- ts(x, start= c(2014, 01), frequency = 36)
#   plot(x, type="l", lty=3)
#     lines(1:length(x), detrend(x), lwd=1.5, col="red") 
#       lines(1:length(x), detrend(x, periodicity = "additive"), 
# 	        lwd=1.5, col="blue") 
# 
# par(mfrow=c(2,2))
#   plot(x, type="l", lty=3, main="raw time-series")
#   plot(detrend(x), type="l", lwd=1.5, col="red", 
#        ylim=range(x), main="Bayesian multi-model")
#   plot(detrend(x, periodicity = "additive"), type="l",  
# 	   lwd=1.5, col="blue", ylim=range(x),
# 	   main="cosine decomposition")
#   plot(x, type="l", lty=3)  
#     lines(1:length(x), detrend(x), lwd=1.5, col="red")  
#       lines(1:length(x), detrend(x, periodicity = "additive"), 
# 	        lwd=1.5, col="blue")  
# 
