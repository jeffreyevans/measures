suppressMessages(
  lapply(c("sf", "spatialEco", "terra", "dplyr", "duckdb", "DBI",  
         "rts", "parallel", "snow", "doSNOW", "xts"), 
		 require, character.only = TRUE))

#***************************
# Set working directory
path = "C:/evans/PFP"
country <- c("brazil", "bhutan", "canada", "colombia", "costa_rica", "peru")[4]

root = file.path(path, country)
setwd(file.path(root, "data"))
  mdl.dir = file.path(root, "model")
  dat.dir = file.path(root, "data")

replace.metric = c(TRUE, FALSE)[1]
mres = c("250m", "300m", "500m")[3]

#*************************************************
# Set intervention date
if(country == "brazil") {
  intervention.date = as.Date("2014-01-10")  
} else if(country == "bhutan") {
  intervention.date = as.Date("2017-01-10")  
} else if(country == "canada") {
  intervention.date = as.Date("2006-01-10")
} else if(country == "colombia") {
  intervention.date = as.Date("2020-01-10")  	
} else if(country == "costa_rica") {
  intervention.date = as.Date("2010-01-10")  
} else if(country == "peru") {
  intervention.date = as.Date("2019-01-10")  
} else {
  warning("Not a valid option")
}  
if(mres == "300m" & intervention.date < as.Date("2015-01-10")) {
  intervention.date = as.Date("2015-01-10")
}


#***********************************************
# LAI causal model response variables
metric.max <- 9
diff.thres = 0.5

m <- rast(paste0("forest_", mres, ".tif"))
  metric <- rast(paste0(toupper(metric.type), "_", mres, "_trend.tif"))
    dates <- as.Date(names(metric))

if(file.exists(file.path(mdl.dir, "lai", "lai_treatment_responses.tif"))){
  r250 <- rast(file.path(mdl.dir, "lai", "lai_treatment_responses.tif"))
  c250 <- rast(file.path(mdl.dir, "lai", "lai_control_responses.tif"))
} else {  
  r250 <- rast(file.path(dat.dir, "forest_250m.tif"))
    r250[] <- NA
      r250 <- rep(r250, 6)
  c250 <- r250
}

lai.con <- file.path(mdl.dir, "lai", "lai_250m_model_data.duckdb")
  con <- dbConnect(duckdb::duckdb(), lai.con, read_only=TRUE)
    db.table <- dbListTables(con)[1]
    dat <- dbGetQuery(con, paste0("SELECT * FROM ", db.table))     
  dbDisconnect(con, shutdown = TRUE)

#if(!any(names(dat) %in% c("gain", "loss"))) {
  current <- dat$current 
  baseline <- dat$pre 
  gain <- current - baseline
    gain[gain == 0] <- 0.01
      gain[gain < 0] <- 0
  loss <- current - baseline
    loss[loss == 0] <- -0.01
      loss[loss > 0] <- 0
        loss <- loss * -1    
  cpct <- (current / baseline) - 1
    lt.idx <- which(cpct < 0)
      for(i in lt.idx) {
        if( abs(current[i] - baseline[i]) <= diff.thres ) {
          cpct[i] <- 0
        } else {
          cpct[i] <- (current[i] / baseline[i]) / metric.max	  
  	  }
      }	
    gt.idx <- which(cpct > 1)
      for(i in gt.idx) {
        if( abs(current[i] - baseline[i]) <= diff.thres ) {
          cpct[i] <- 0
        } else {
          cpct[i] <- (current[i] - baseline[i]) / metric.max	  
  	    }
      }
  summary(gain); summary(loss); summary(cpct)
  dat$gain <- gain 
  dat$loss <- loss
  dat$cpct <- cpct

#********************************************************
# Porportion timeseires is over/under baseline
metric <- rast("lai_250m__trend.tif"))
  dates <- as.Date(names(metric))
    pre.idx <- which(dates < intervention.date)   
    post.idx <- which(dates >= intervention.date)   

gain.pct <- vector()
loss.pct <- vector()
  for(i in 1:nrow(dat)) {
    cat("iteration", i, "of", nrow(dat), "\n")
    x <- as.numeric(metric[dat[i,]$cell])[post.idx]
	loss.pct[i] <- sum(x <= dat[i,]$pre) / length(post.idx)
	gain.pct[i] <- sum(x > dat[i,]$pre) / length(post.idx) 
  }

# library(future.apply)
# plan(multisession)
# 
# gain <- future_lapply(1:nrow(dat), FUN = \(o) {
#     return(sum(as.numeric(metric[dat[o,]$cell])[post.idx] > dat[o,]$pre) / length(post.idx))
#   })

  
  # Write updat  gain.loss(x, p = dat[i,]$pre)
  con <- dbConnect(duckdb::duckdb(), lai.con, read_only=FALSE)
    db.table <- dbListTables(con)[1]
    duckdb::dbWriteTable(con, paste0(db.table), dat, overwrite = TRUE)
  dbDisconnect(con, shutdown = TRUE)   	   
#}

# split treatment and control pixels
nidx <- which(dat$NAME %in% "control")
  if(length(nidx) > 0) {
    ctl <- dat[nidx,]
    dat <- dat[-nidx,]
  }

# Write rasters
r250[dat$cell] <- dat[,c("pre", "current", "tsa", "cpct", "gain", "loss")]
  names(r250) <- c("pre", "current", "tsa", "cpct", "gain", "loss") 
    if(file.exists(file.path(mdl.dir, "lai", "lai_treatment_responses.tif"))){
      stars::write_stars(stars::st_as_stars(r250), dsn = file.path(mdl.dir, "lai", 
	                     "lai_treatment_responses.tif"), type = "Float32", 
						 update = TRUE)
    } else {
      writeRaster(r250, file.path(mdl.dir, "lai", "lai_treatment_responses.tif"),
                  overwrite=TRUE) 
    }

c250[ctl$cell] <- ctl[,c("pre", "current", "tsa", "cpct", "gain", "loss")]
  names(c250) <- c("pre", "current", "tsa", "cpct", "gain", "loss") 
    if(file.exists(file.path(mdl.dir, "lai", "lai_control_responses.tif"))){
      stars::write_stars(stars::st_as_stars(c250), dsn = file.path(mdl.dir, "lai", 
	                     "lai_control_responses.tif"), type = "Float32", 
						 update = TRUE)
    } else {
      writeRaster(c250, file.path(mdl.dir, "lai", "lai_control_responses.tif"),
                  overwrite=TRUE) 
    }

remove(dat, current, baseline, gain, loss, cpct, r250, c250)

#***********************************************
# fCOV causal model response variables
metric.max <- 1
diff.thres = 0.01

if(file.exists(file.path(mdl.dir, "fcov", "fcov_treatment_responses.tif"))){
  r500 <- rast(file.path(mdl.dir, "fcov", "fcov_treatment_responses.tif"))
  c500 <- rast(file.path(mdl.dir, "fcov", "fcov_control_responses.tif"))
} else {  
  r500 <- rast(file.path(dat.dir, "forest_500m.tif"))
    r500[] <- NA
      r500 <- rep(r500, 6)
  c500 <- r500
}

fcov.con <- file.path(mdl.dir, "fcov", "fcov_500m_model_data.duckdb")
  con <- dbConnect(duckdb::duckdb(), fcov.con, read_only=TRUE)
    db.table <- dbListTables(con)  # list tables in Db
    fdat <- dbGetQuery(con, paste0("SELECT * FROM ", db.table))     
  dbDisconnect(con, shutdown = TRUE)

#if(!any(names(fdat) %in% c("gain", "loss", "cpct"))) {
  current <- fdat$current 
  baseline <- fdat$pre 
  gain <- current - baseline
    gain[gain == 0] <- 0.01
      gain[gain < 0] <- 0
  loss <- current - baseline
    loss[loss == 0] <- -0.01
      loss[loss > 0] <- 0
        loss <- loss * -1    
  cpct <- (current / baseline) - 1
    lt.idx <- which(cpct < 0)
      for(i in lt.idx) {
        if( abs(current[i] - baseline[i]) <= diff.thres ) {
          cpct[i] <- 0
        } else {
          cpct[i] <- (current[i] / baseline[i]) / metric.max	  
  	  }
      }	
    gt.idx <- which(cpct > 1)
      for(i in gt.idx) {
        if( abs(current[i] - baseline[i]) <= diff.thres ) {
          cpct[i] <- 0
        } else {
          cpct[i] <- (current[i] - baseline[i]) / metric.max	  
  	    }
      }
  summary(gain); summary(loss); summary(cpct)	
  fdat$gain <- gain 
  fdat$loss <- loss
  fdat$cpct <- cpct 	

  # Write updated database table
  con <- dbConnect(duckdb::duckdb(), fcov.con, read_only=FALSE)
    db.table <- dbListTables(con)[1]
    duckdb::dbWriteTable(con, paste0(db.table), fdat, overwrite = TRUE)
  dbDisconnect(con, shutdown = TRUE)   	
#}

# split treatment and control pixels
nidx <- which(fdat$NAME %in% "control")
  if(length(nidx) > 0) {
    fctl <- fdat[nidx,]
    fdat <- fdat[-nidx,]
  }

# Write rasters
r500[fdat$cell] <- fdat[,c("pre", "current", "tsa", "cpct", "gain", "loss")]
  names(r500) <- c("pre", "current", "tsa", "cpct", "gain", "loss")
    if(file.exists(file.path(mdl.dir, "fcov", "fcov_treatment_responses.tif"))){
      stars::write_stars(st_as_stars(r500), dsn = file.path(mdl.dir, "fcov", "fcov_treatment_responses.tif"), 
                         type = "Float32", update = TRUE)
    } else {
      writeRaster(r500, file.path(mdl.dir, "fcov", "fcov_treatment_responses.tif"),
                  overwrite=TRUE) 
    }

c500[fctl$cell] <- fctl[,c("pre", "current", "tsa", "cpct", "gain", "loss")]
  names(c500) <- c("pre", "current", "tsa", "cpct", "gain", "loss")
    if(file.exists(file.path(mdl.dir, "fcov", "fcov_control_responses.tif"))){
      stars::write_stars(st_as_stars(c500), dsn = file.path(mdl.dir, "fcov", "fcov_control_responses.tif"), 
                         type = "Float32", update = TRUE)
    } else {
      writeRaster(c500, file.path(mdl.dir, "fcov", "fcov_control_responses.tif"),
                  overwrite=TRUE) 
    }

remove(fdat, current, baseline, gain, loss, cpct, r500, c500)

# 
# 
# # Plot distributions
# treat.cpct.den <- density(as.numeric(na.omit(m500[["cpct"]][])[,1]), bw="SJ")
# ctl.cpct.den <- density(as.numeric(na.omit(c500[["cpct"]][])[,1]), bw="SJ")
# 
# plot(treat.cpct.den, type="n")
#   polygon(treat.cpct.den, col=rgb(1, 0, 0, 0.5))
#   polygon(ctl.cpct.den, col=rgb(0, 0, 1, 0.5))
# 
# #*************************************************
# # Summarize effect sizes by ELU/WTE
# lai.es <- rast("lai_effect_sizes.tif")
#   lai.mask <- sum(lai.es, na.rm=TRUE)
# fcov.es <- rast("fcov_effect_sizes.tif")
#   fcov.mask <- sum(fcov.es, na.rm=TRUE)
# 
# wte250 <- mask(rast(file.path(dat.dir, "elu_250m.tif")), lai.mask)
# wte500 <- mask(rast(file.path(dat.dir, "elu_500m.tif")), fcov.mask)
# 
# lai.sum <- zonal(lai.es, wte250, summary, na.rm=TRUE)
# f250 <- freq(wte250)
#   f250 <- tapply(f250$count, f250$value, sum)
#     lai.sum$counts[which(lai.sum$elu %in% names(f250))] <- as.numeric(f250)
# write.csv(lai.sum, file.path(getwd(), "tables", "lai_wte_summary.csv"))  
#   
# fcov.sum <- zonal(fcov.es, wte500, summary, na.rm=TRUE)
#   f500 <- freq(wte500)
#   f500 <- tapply(f500$count, f500$value, sum)
#     fcov.sum$counts[which(fcov.sum$elu %in% names(f500))] <- as.numeric(f500)
# write.csv(fcov.sum, file.path(getwd(), "tables", "fcov_wte_summary.csv"))  
# 