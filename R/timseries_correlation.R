library(terra)
library(spatialEco)
library(sf)
library(future.apply)
library(bfast)
library(imputeTS)

root = "C:/evans/Colombia/waterfunds"
  dat.dir <- file.path(root, "data")
setwd(root)

source(file.path(root, "code", "raster.kendall.R"))

mdl.dir <- file.path(root, "model", "correlation")
  metrics <- c("lai", "fcov")[2]

#***************************************
# read treatment units 
# st_layers(file.path(dat.dir, "CuencaVerde.gpkg"))
treat <- st_read(file.path(dat.dir, "CuencaVerde.gpkg"), "treatments")

#***************************************
# Calculate Mann-Kendall TAU and Sen Slope
#
# read raster time-series and parse dates
#dry.season <- which(months(dates) %in% unique(months(dates))[1:4]) 
#wet.season <- which(months(dates) %in% unique(months(dates))[5:12]) 
metric <- rast(file.path(dat.dir, paste0(toupper(metrics), "_300m_trend", ".tif")))
  d <- unlist(lapply(strsplit(names(metric), "_"), function(x) x[2]))
  dates <- as.Date(unlist(lapply(d, function(j) {
      paste( substr(j,1,4), substr(j,5,6), substr(j,7,8), sep="-" )
    })))

metric.tau <- raster.kendall(metric)
  tau.pos <- ifel(metric.tau$slope > 0, metric.tau$tau, NA)
    names(tau.pos) <- "tau_positive" 
  tau.neg <- ifel(metric.tau$slope < 0, metric.tau$tau, NA)
    names(tau.neg) <- "tau_negative" 
 metric.tau <- c(metric.tau, tau.pos, tau.neg)  	
  terra::writeRaster(metric.tau, file.path(mdl.dir, 
                     paste0(tolower(metrics), "_tau", ".tif")),
                     overwrite = TRUE)
			  
plot(metric.tau[["tau"]])
  plot(st_geometry(treat),add=TRUE)

# *************************************
# Autocorrelation corrected Kendall Tau 
#   with detrended timeseires
mk.tau <- function(x, min.obs = 10, ...) {
  if(length(x[!is.na(x)]) < min.obs) {
    return( rep(NA,5) )
  } else {
   return( round(spatialEco::kendall(x[!is.na(x)], ...),5) )
  }
}
# metric.tau <- terra::app(metric, mk.tau, cores = 4) 
tau.dat <- as.data.frame(metric, cells=TRUE, na.rm=NA)
  cell.ids <- tau.dat[,1] 
    tau.dat <- tau.dat[,-1]

system.time({ 
  zyp::zyp.trend.vector(as.numeric(tau.dat[10,]), method = "zhang")
})

system.time({ 
  mk.tau(as.numeric(tau.dat[10,]), method = "zhang")
})

system.time({ 
  tau <- apply(tau.dat, 1, mk.tau, method="zhang") 
})
  # tau <- as.data.frame(t(tau))

# tau.dat <- split(tau.dat, seq(nrow(tau.dat)))
# system.time({ 
#   tau <- lapply(tau.dat, function(x) mk.tau(as.numeric(x[1,]), method="none")) 
# })
# tau <- do.call(rbind, tau)

metric.tau <- metric[[1]]
  metric.tau[] <- rep(NA, ncell(metric.tau))
    metric.tau <- rep(metric.tau, 6)
      names(metric.tau) <- names(tau)

metric.tau[cell.ids] <- tau
  tau.pos <- ifel(metric.tau$slope > 0, metric.tau$tau, NA)
    names(tau.pos) <- "tau_positive" 
  tau.neg <- ifel(metric.tau$slope < 0, metric.tau$tau, NA)
    names(tau.neg) <- "tau_negative"  
metric.tau <- c(metric.tau, tau.pos, tau.neg)

terra::writeRaster(metric.tau, file.path(mdl.dir, 
                    paste0(tolower(metrics), "_tau_autocor", ".tif")),
                    overwrite = TRUE)
