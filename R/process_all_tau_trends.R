
#***************************************
#***************************************
# process countries
for(idx in 1:3) { 
  for(m in c("lai", "fcov")) { 

    if(m == "lai"){
      mres = "250m" 
    } else if(m == "fcov") {		
      mres = "500m" 
    }
 
    cty <- c("bhutan", "costa_rica", "canada")[idx]
    root = file.path("Z:/PFP", cty)
    setwd(root)
      dat.dir = file.path(root, "data")
      tau.dir = file.path(root, "model", m, "tau")

    cat("****", "\n")
    cat("**** Processing", m, "trends for", cty, "\n")  
    cat("****", "\n")
    flush.console(); Sys.sleep(0.01)
	
    #***************************************
    #***************************************
    # timeseries
    metric <- rast(file.path(dat.dir, paste0(m, "_", mres, "_trend.tif")))
      dates <- as.Date(names(metric))
	  #dates <- as.Date(names(rast("Z:/PFP/colombia/data/lai_250m_trend.tif")))
      
    # if n > max.size, create index for iteration
    n <- nrow(metric) * ncol(metric)
      if(n < max.size) {
        a <- list(range(1:n))
      } else { 
        a <- lapply(split(1:n, ceiling(seq_along(1:n) / s)), range)
      }

    # column indices for pre/post dates
    date.break = c("2017-01-01", "2010-01-01", "2006-01-01", "2020-01-01", "2019-01-01", "2014-01-01")[idx]
      pre.idx <- which(dates < date.break) 
      post.idx <- which(dates >= date.break) 
    
    # create empty output rasters
    forest <- rast(file.path(dat.dir, paste0("forest_", mres, ".tif")))
      forest[forest == 0] <- NA
    f <- rast(file.path(dat.dir, paste0("forest_", mres, ".tif")))
      f[] <- NA

    out.pre <- file.path(tau.dir, paste0("full_", m, "_pre_tau_" , mres, ".tif"))    # output pre intevention tau raster results
    out.post <- file.path(tau.dir, paste0("full_", m, "_post_tau_" , mres, ".tif"))  # output pre intevention tau raster results
     
    pre.tau <- rep(f, L)
      names(pre.tau) <- tau.names
        pre.tau <- writeRaster(pre.tau, out.pre, datatype="FLT4S", overwrite=TRUE)
    post.tau <- rep(f, L)
      names(post.tau) <- tau.names
        post.tau <- writeRaster(post.tau, out.post, datatype="FLT4S", overwrite=TRUE)
    	
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
        cells <- start.row:end.row 
          tau.dat <- metric[cells]
            na.idx <- which(is.na(tau.dat[,100])) 
              if(length(na.idx) > 0) tau.dat <- tau.dat[-na.idx,]
                if(nrow(tau.dat) < 1) next
    	  cells <- cells[-na.idx]
    	  
    	# TAU pre-intervention  	 
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
          pre.tau[cells] <- tau.pre
            pre.tau <- mask(pre.tau, forest)	  
    	      stars::write_stars(st_as_stars(pre.tau), dsn = out.pre,
                                 type = "Float32", update = TRUE)
    	remove(tau.pre)
    	  
    	# TAU post-intervention  	  
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
          post.tau[cells] <- tau.post
    	    post.tau <- mask(post.tau, forest)
    	      stars::write_stars(st_as_stars(post.tau), dsn = out.post,
                                 type = "Float32", update = TRUE)
        remove(tau.post, tau.dat, start.row, end.row)	
      gc()  
    }
  }
}

## Clean up all raster temp files
tmpFiles(current=TRUE, orphan=TRUE, old=TRUE, remove=TRUE)

