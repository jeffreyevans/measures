# Performs a Kendall trend analysis for lai or fcov timeseries
#   subsets pre/post intervention following;
# PFP                   Start   Year range
# Great Bear	        2006	2005-2021
# Forever Costa Rica   	2010	2009-2021
# ARPA					2014	2013-2021
# Bhutan for Life		2017	2016-2021
# Peru					2019	2018-2021
# Colombia 			**	2022	2020-2021
#
# Allows for serial autocorrelation correction (prewhitening) using the Yue & Pilon (2002) 
# method or uses a standard Kendall test. Temporal data can be colloposed into quarters or years
# for more efficent processing of long timeseries.  
# 
# Yue, S., P. Pilon, B. Phinney and G. Cavadias, 2002. The influence of autocorrelation on the ability 
#   to detect trend in hydrological series. Hydrological Processes, 16: 1807-1829. 
#
suppressMessages(
  lapply(c("sf", "spatialEco", "terra", "dplyr", "duckdb", "DBI",  
         "rts", "fst", "parallel", "snow", "doSNOW", "xts"), 
		 require, character.only = TRUE))

#***************************************
#***************************************
# set up environment
m = c("lai", "fcov")[2]
mres = c("250m", "300m", "500m")[3]
  if(mres == "250m") { 
    out.res = 250
    m = "lai"
  } else if(mres == "500m") { 
    out.res = 500
    m = "fcov"		
  } else { 
    out.res = 300 
  }

idx = 4  # indicates country 
cty <- c("bhutan", "costa_rica", "canada", "colombia", "peru", "brazil")[idx]
root = file.path("C:/evans/PFP", cty)

setwd(root)
  dat.dir = file.path(root, "data")
  mdl.dir = file.path(root, "model", paste0("model",mres))
  tau.dir = file.path(mdl.dir, "tau")
  db.dir = file.path(mdl.dir, "database")
dir.create(tau.dir, showWarnings = FALSE)
  
max.size = 250000                             # maximum table size before iteration
s = 250000                                    # subset sizes in iteration
autocorrelation = c(TRUE, FALSE)[1]           # correct for serial autocorrelation (Yuepilon correction)
np = future::availableCores()-1               # number of processing cores
agg.data = c("none", "quarters", "years")[2]  # aggregrate data to quarters or years

# database connections and raster output
ts.db <- file.path(db.dir, paste0(m, "_",mres, "_timeseries_data.duckdb"))           # timeseries data connection
dat.db <- file.path(mdl.dir, paste0(m, "_",mres, "_model_data.duckdb"))           # timeseries data connection
tau.db <- file.path(tau.dir, paste0(m, "_tau_" , mres, ".duckdb"))                   # output tau database results
out.pre <- file.path(tau.dir, paste0(m, "_pre_tau_" , mres, ".tif"))                 # output pre intevention tau raster results
out.post <- file.path(tau.dir, paste0(m, "_post_tau_" , mres, ".tif"))               # output pre intevention tau raster results

# Date ranges
date.break = c("2017-01-01", "2010-01-01", "2006-01-01", "2020-01-01", "2019-01-01", "2014-01-01")[idx]
ts.date <- list()
  ts.date[[1]] <- c(2017, 1, 1)
    ts.date[[2]] <- c(2010, 1, 1)
      ts.date[[3]] <- c(2006, 1, 1)
    ts.date[[4]] <- c(2020, 1, 1)
  ts.date[[5]] <- c(2019, 1, 1)
ts.date[[6]] <- c(2014, 1, 1)

if(autocorrelation) {
  test = "yuepilon"
  L=7 
  tau.names <- c("tau", "trend", "slope", "p_value", "intercept", "LCL", "UCL") 

} else { 
  test = "mk.cor"
  L=3
  tau.names <- c("tau", "z_value", "p_value")  
}

#***************************************
#***************************************
# functions
mk.tau <- function(x, min.obs = 5, offset = 2, d = NULL,  
                   ts.period = c("none", "years", "quarters", "months", "weeks"), 
                   type = c("mk.cor", "kendall", "lm", "yuepilon", "zhang")) {
	knames <- c("tau", "trend", "slope", "intercept", "p_value", "LCL", "UCL", "z_value")
	ts.period = ts.period[1]
	  type = type[1]
    x <- as.numeric(x)
	if(ts.period != "none") {
	  if(is.null(d) | length(d) != length(x)) { 
	    warning("Cannot aggregrate by specified timeperiod")
      } else {
        n = length(x)	  
        x <- xts::xts(x, order.by = d)
        x <- as.data.frame(xts::to.period(x, period = ts.period))
	    d <- as.Date(rownames(x))
	    x <- apply(x, MARGIN=1, FUN=mean)
	    message("Aggregrated n=", n, " into n=", length(x), " ", ts.period, " observations") 	
	  }
	}
  if(type %in% c("zhang", "yuepilon")) { 
    n = 7 
  } else if(type == "mk.cor") {
    n = 3   
  } else if(type %in% c("lm", "kendall")) {
    n = 6
  }
  if(length(x[!is.na(x)]) < min.obs) {
    tau <- rep(NA,n) 
    } else { 
      if(type == "Zhang") {
		kcol <- knames[c(1,2,3,5,4,6,7)]
  	      tau <- round(zyp::zyp.zhang(x[!is.na(x)]),5)[c(5,3,2,11,6,1,4)]
  		    names(tau) <- kcol 
      }  else if(type == "yuepilon") {
		kcol <- knames[c(1,2,3,5,4,6,7)]
  	      tau <- round(zyp::zyp.yuepilon(x[!is.na(x)]),5)[c(5,3,2,11,6,1,4)]
  		    names(tau) <- kcol 
      }  else if(type == "mk.cor") {
        kcol <- knames[c(1,8,5)]		
		  tau <- round(as.data.frame(stats::cor.test(1:length(x[!is.na(x)]), x[!is.na(x)], 
		               method = "kendall")[c(4,1,3)]),6)
          names(tau) <- kcol
      }  else if(type == "kendall") {
	    kcol <- knames[c(1,3,4,5,6,7)]	
          tau <- round(spatialEco::kendall(x[!is.na(x)]),5)[c(2,1,3,4,5,6)] 
            names(tau) <- kcol 
      }  else if(type == "lm") {
		kcol <- knames[c(1,3,4,5,6,7)]	 
		  xts <- as.numeric(x[-c(1:offset)])
		  yts <- as.numeric(x[-c((length(x)-(offset-1)):length(x))]) 
            p <- stats::predict(stats::lm(xts ~ yts))		
              tau <- round(spatialEco::kendall(p),5)[c(2,1,3,4,5,6)]
			    names(tau) <- kcol 	 
	  }
  }
  return( tau )
}

#***************************************
#***************************************
# timeseries database
con <- dbConnect(duckdb::duckdb(dbdir = ts.db), read_only = FALSE)
  db.table <- dbListTables(con)[1]
  n <- dbGetQuery(con, paste0("SELECT COUNT(*) FROM ", db.table))[,1]
  ts.names <- dbGetQuery(con, paste0("DESCRIBE ", db.table))[,"column_name"] 
dbDisconnect(con, shutdown = TRUE)

if(n < max.size) {
  a <- list(range(1:n))      	
} else { 
  a <- lapply(split(1:n, ceiling(seq_along(1:n) / s)), range)
}

rm.cols <- which(ts.names %in% c("cell","elu"))

# Create dates and indices
d <- gsub("D", "", ts.names[-c(1:2)])
dates <- as.Date(unlist(lapply(d, function(j) {
            paste( substr(j,1,4), substr(j,5,6), substr(j,7,8), sep="-" )
          })))
pre.idx <- which(dates < date.break)
post.idx <- which(dates >= date.break)

# create empty output rasters
f <- rast(file.path(dat.dir, paste0("forest_", mres, ".tif")))
  f[] <- NA
pre.tau <- rep(f, L)
  names(pre.tau) <- tau.names
post.tau <- rep(f, L)
  names(post.tau) <- tau.names

#*****************************************************************
#*****************************************************************
# pre and post intervention trend based on autocorrelation corrected 
# decorrelated timeseires using ARIMA
#*****************************************************************
#*****************************************************************
for(i in 1:length(a)) {
  cat("iteration", i, "of", length(a), "\n") 
  start.row = a[[i]][1] 
  end.row = a[[i]][2]
  con <- dbConnect(duckdb::duckdb(dbdir = ts.db), read_only = FALSE)
    db.table <- dbListTables(con) 
    tau.dat <- dbGetQuery(con, paste0("SELECT * FROM ", db.table, " WHERE rowid BETWEEN ", start.row, " AND ", end.row ))       
  dbDisconnect(con, shutdown = TRUE)   
    cells <- tau.dat$cell
	  tau.dat <- tau.dat[,-rm.cols]  
  con <- dbConnect(duckdb::duckdb(dbdir = dat.db), read_only = FALSE)
    db.table <- dbListTables(con) 
    dat <- dbGetQuery(con, paste0("SELECT * FROM ", db.table, " WHERE rowid BETWEEN ", start.row, " AND ", end.row ))       
  dbDisconnect(con, shutdown = TRUE) 
    dat <- dat[,c(1:14)]
      ctl.idx <- which(dat$NAME == "control")
    cat(cty, m, mres, "pre-intervention data block has", format(nrow(tau.dat), big.mark=",", 
        scientific=FALSE, trim=TRUE), "obs and", length(pre.idx), "timeseries", "\n")
        flush.console(); Sys.sleep(0.01)
      system.time({		
        cl <- makeCluster(detectCores()-1, type = "SOCK", outfile="")
    	    clusterEvalQ(cl, {library(stats)})
              clusterExport(cl, c("dates", "mk.tau", "pre.idx", "test", "agg.data"), 
			                envir=environment())
              registerDoSNOW(cl)
                tau.pre <- parApply(cl, tau.dat[,pre.idx], MARGIN=1, FUN = \(o) {
      		                        as.data.frame(t(mk.tau(o, d = dates[pre.idx], 
							        type = test, ts.period = agg.data))) } )	 
          stopCluster(cl)
        registerDoSEQ()
    })
    tau.pre <- do.call(rbind, tau.pre)
      pre.tau[dat$cell[-ctl.idx]] <- tau.pre[-ctl.idx,] 
	    tau.pre <- data.frame(period = paste(c("pre", as.character(range(dates[pre.idx]))), collapse = " "),
                              dat, tau.pre)	
	con <- dbConnect(duckdb::duckdb(dbdir = tau.db), read_only = FALSE)
      if(!DBI::dbExistsTable(con, "tau")) {
        duckdb::dbWriteTable(con, "tau", tau.pre)
        dbDisconnect(con, shutdown = TRUE)
      } else {
        duckdb::dbAppendTable(con, "tau", tau.pre)
        dbDisconnect(con, shutdown = TRUE)
      }
	  remove(tau.pre) 
    cat(cty, m, mres, "post-intervention data block has", format(nrow(tau.dat), big.mark=",", 
        scientific=FALSE, trim=TRUE), "obs and", length(post.idx), "timeseries", "\n")
        flush.console(); Sys.sleep(0.01)
      system.time({		
        cl <- makeCluster(detectCores()-1, type = "SOCK", outfile="")
    	    clusterEvalQ(cl, {library(stats)})
              clusterExport(cl, c("dates", "mk.tau", "post.idx", "test", "agg.data"), 
			                envir=environment())
              registerDoSNOW(cl)
                tau.post <- parApply(cl, tau.dat[,post.idx], MARGIN=1, FUN = \(o) {
       		                         as.data.frame(t(mk.tau(o, d = dates[post.idx], 
							         type = test, ts.period = agg.data))) } )	 
          stopCluster(cl)
        registerDoSEQ()
    })
    tau.post <- do.call(rbind, tau.post)
      post.tau[dat$cell[-ctl.idx]] <- tau.post[-ctl.idx,] 
	    tau.post <- data.frame(period = paste(c("post", as.character(range(dates[post.idx]))), collapse = " "),
                     dat, tau.post)
	con <- dbConnect(duckdb::duckdb(dbdir = tau.db), read_only = FALSE)
      if(!DBI::dbExistsTable(con, "tau")) {
        duckdb::dbWriteTable(con, "tau", tau.post)
        dbDisconnect(con, shutdown = TRUE)
      } else {
        duckdb::dbAppendTable(con, "tau", tau.post)
        dbDisconnect(con, shutdown = TRUE)
      }
    remove(tau.post, tau.dat, start.row, end.row, ctl.idx)	
  gc()  
}

writeRaster(pre.tau, out.pre, overwrite=TRUE, datatype="FLT4S")
writeRaster(post.tau, out.post, overwrite=TRUE, datatype="FLT4S")

## Clean up all raster temp files
tmpFiles(current=TRUE, orphan=TRUE, old=TRUE, remove=TRUE)
