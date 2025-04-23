# process time series filling NA's, optionally smoothing 
# and decomposing periodicity to derive a trend. This uses
# a linear trend method or multi-model Bayesian averaging 
# via the BEAST method
#
# in:  metric_res_raw.tif
# out: metric_res_trend.tif 
#
suppressMessages(
  lapply(c("terra", "Rbeast", "xts", "parallel", "snow", "doSNOW", 
           "duckdb", "DBI"), require, character.only = TRUE))

idx = 6  # indicates country
country <- c("bhutan", "costa_rica", "canada", "colombia", "peru", "brazil")[idx]

metric.type <- c("lai", "fcov")[2]
mres = c("250m", "300m", "500m")[3]
  if(mres == "250m") metric.type = "lai"
  if(mres == "500m") metric.type = "fcov"

out.file <- paste0(toupper(metric.type), "_", mres, "_trend.tif") 
in.file <- paste0(metric.type, "_", mres, "_raw.tif")  

root = file.path("Z:/PFP", country)
  dat.dir = file.path(root, "data")
  mdl.dir = file.path(root, "model", paste0("model",mres))
setwd(dat.dir)

max.size = 500000   # maximum matrix size before iteration is implemented
b = 500000          # iteration bin size

build.db = c(TRUE,FALSE)[1] 
collapse.ts = c(TRUE,FALSE)[1]
ts.period = c("none", "quarters", "months")[3]
np = future::availableCores()-1               # number of processing cores
ts.db <- file.path(mdl.dir, paste0(metric.type, "_",mres, "_timeseries_data.duckdb"))    # timeseries data

#*************************************************
# Read metric and forest mask
metric <- rast(file.path(dat.dir, in.file))
  if(mres == "250m" | mres == "500m") { 
   dates <- as.Date(names(metric))
  } else {
    d <- unlist(lapply(strsplit(names(metric), "_"), function(x) x[2]))
    dates <- as.Date(unlist(lapply(d, function(j) {
        paste( substr(j,1,4), substr(j,5,6), substr(j,7,8), sep="-" )
      })))
  }

m <- rast(file.path(dat.dir, paste0("forest_", mres, ".tif")))
 m[m == 0] <- NA

# forest cell indices
idx <- as.data.frame(m, cells = TRUE, na.rm = TRUE)$cell

# Create cell vectors for iteration
n = length(idx)
if(n < max.size) {
  s <- list(idx)
} else { 
  s <- split(idx, ceiling(seq_along(idx) / b))
}

# Create empty raster timeseries
r.trend <- metric

if(collapse.ts) {
  x <- as.numeric(metric[idx[1]])
  x <- as.data.frame(xts::to.period(xts::xts(x, order.by = dates), period = ts.period))
  sub.dates <- as.Date(rownames(x))
  r.trend <- r.trend[[1:length(sub.dates)]]
    names(r.trend) <- sub.dates
	d = sub.dates
  remove(x)
} else {
  d <- dates
}

#*************************************************
# calculate multi-threaded BEAST seasional detrend model

ts.collapse <- function(x, p = ts.period, d = dates) {
  x <- xts::to.period(xts::xts(x, order.by = d), period = p)
  as.numeric(apply(x, MARGIN=1, FUN=mean))
}  
	   
# BEAST model parameters 
extra.args <- list()
  extra.args$quiet = TRUE
  extra.args$printOptions = FALSE
  extra.args$computeTrendSlope = TRUE
  extra.args$printProgressBar = TRUE
parms = list()        
  parms$whichDimIsTime = 2
  parms$maxMissingRateAllowed = 0.50 
  parms$period = 1.0
  parms$start=dates[1]
  #parms$time = dates
  if(mres == "250m" | mres == "500m") {
    parms$startTime = dates[1]
    parms$deltaTime = 1/46 
  } else if(mres == "300m") {
	parms$startTime = dates[1]
    parms$deltaTime = 1/20
 }

# iteration for cell indices
  for(j in 1:length(s)) {
    i = s[[j]]
	cat(j, "of", length(s), "iterations", "n =", length(i), "\n")
    tsm <- metric[i]
      if(collapse.ts) {
	    cat("Collapsing timeseries to", ts.period, "\n") 
        cl <- makeCluster(detectCores()-1, type = "SOCK", outfile="")
          clusterEvalQ(cl, {library(xts)})
            clusterExport(cl, c("ts.collapse", "dates", "ts.period"), 
  			            envir=environment())
                registerDoSNOW(cl)
                  tsm <- parApply(cl, tsm, MARGIN=1, FUN = \(o) { ts.collapse(o) } )
            stopCluster(cl)
          registerDoSEQ()
  	    tsm <- t(tsm)
		   parms$startTime = sub.dates[1]
           parms$deltaTime = ifelse(ts.period == "quarters", 1/4, 
                               ifelse(ts.period == "months", 1/12, NA))
	  gc()
      }
    cat("Using Bayesian (BEAST) method to detrend", format(nrow(tsm), big.mark=",", 
      scientific=FALSE, trim=TRUE), "obs and", ncol(tsm), "timeseries", "\n")
        flush.console(); Sys.sleep(0.01) 
      o <- beast123(as.matrix(tsm), metadata = parms, season  = "harmonic", 
    		        method  = "bayes", extra = extra.args)
      mm <- c(min(o$trend$Y[!is.nan(o$trend$Y)], na.rm=TRUE), 
    		    max(o$trend$Y[!is.nan(o$trend$Y)], rm=TRUE))  
    	if(mm[1] < 0){ o$trend$Y[o$trend$Y < 0] <- 0 } 
        if(metric.type == "fcov"& mm[2] > 1) {
    	    o$trend$Y[o$trend$Y > 1] <- 1      
  	    }	  
        if(metric.type == "lai" & mm[2] > 9) {
    	    o$trend$Y[o$trend$Y > 9] <- 9 
    	}
		o$trend$Y <- round(o$trend$Y, 4)
      r.trend[i] <- as.data.frame(o$trend$Y)
    if(write.db) {
	  tsm <- data.frame(cell=i,o$trend$Y) 
	    dnames <- paste0("D", gsub("-", "", as.character(d)))
	      names(tsm)[-1] <- dnames
      cat("Opening database connection and writing writing timeseries", "\n")
        con <- dbConnect(duckdb::duckdb(dbdir = ts.db), read_only = FALSE)
          if(!DBI::dbExistsTable(con, metric.type)) {
            duckdb::dbWriteTable(con, metric.type, tsm, append = FALSE)
            dbDisconnect(con, shutdown = TRUE)
          } else {
            duckdb::dbAppendTable(con, metric.type, tsm)
            dbDisconnect(con, shutdown = TRUE)
          }	  
	}   
    remove(tsm, i, o, mm)
  gc() 
  }

terra::writeRaster(r.trend, out.file, overwrite=TRUE, datatype="FLT4S")

terra::tmpFiles(current=TRUE, remove=TRUE)
