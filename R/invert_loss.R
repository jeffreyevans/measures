suppressMessages(
  lapply(c("sf", "spatialEco", "terra",  "dplyr", "duckdb", 
           "DBI", "grf", "ggplot2", "data.table"), 
			require, character.only = TRUE))

clrs <- c("darkred","darkorange1","azure3","aquamarine3","darkgreen")
bks <- c( seq(-1,-0.005, 0.10), c(-0.01, 0.01), seq(0.10, 1, 0.1)) 
hte.mcls <- colorRampPalette(clrs)(21)     

#***************************
# Set working directory
path = "C:/evans/PFP"
country <- c("brazil", "bhutan", "canada", "colombia", "costa_rica", "peru")[5]

root = file.path(path, country)
setwd(file.path(root, "data"))
  mdl.dir = file.path(root, "model")
  dat.dir = file.path(root, "data")

#******************************************************
# Correct LAI effect sizes for zero non gains or losses
#******************************************************
lai.names <- c(list.files(file.path(mdl.dir, "lai", "gain"), "tif$", full.names=TRUE),
               list.files(file.path(mdl.dir, "lai", "loss"), "tif$", full.names=TRUE))

lai.con <- file.path(mdl.dir, "lai", "lai_250m_model_data.duckdb")
  con <- dbConnect(duckdb::duckdb(), lai.con, read_only=TRUE)
    db.table <- dbListTables(con)[1]
    dat <- dbGetQuery(con, paste0("SELECT * FROM ", db.table))     
  dbDisconnect(con, shutdown = TRUE)

nidx <- which(dat$NAME %in% "control")
  if(length(nidx) > 0) { dat <- dat[-nidx,] }

#***************************
# Correct loss with no changes
neg.idx <- which(dat$diff >= -0.5 & dat$diff <= 0.5)
  summary(dat$diff[neg.idx])
neg.cells <- dat$cell[neg.idx]

lai.loss <- rast(grep("loss", lai.names, value = TRUE))
  lai.loss[neg.cells] <- 0  

# write raster
loss.out <- file.path(dirname(grep("loss", lai.names, value = TRUE)),
      paste0(rm.ext(basename(grep("loss", lai.names, value = TRUE))), 
	  "_corrected.tif"))
writeRaster(lai.loss, loss,out, overwrite = TRUE) 

plot(lai.loss[[1]], maxcell=5000000, smooth = TRUE, type="interval",  
     breaks = bks, col = hte.mcls)
       title("lai loss effect size")

#***************************
# Correct gain with no changes
pos.idx <- which(dat$diff >= -0.5 & dat$diff <= 0.5)
  summary(dat$diff[pos.idx])
pos.cells <- dat$cell[pos.idx]

lai.gain <- rast(grep("gain", lai.names, value = TRUE))
  lai.gain[pos.cells] <- 0  

# write raster
gain.out <- file.path(dirname(grep("gain", lai.names, value = TRUE)),
      paste0(rm.ext(basename(grep("gain", lai.names, value = TRUE))), 
	  "_corrected.tif"))
writeRaster(lai.gain, gain.out, overwrite = TRUE) 

plot(lai.gain[[1]], maxcell=5000000, smooth = TRUE, type="interval",  
     breaks = bks, col = hte.mcls)
       title("lai gain effect size")

#******************************************************
# Correct LAI effect sizes for zero non gains or losses
#******************************************************
lai.names <- c(list.files(file.path(mdl.dir, "lai", "gain"), "tif$", full.names=TRUE),
               list.files(file.path(mdl.dir, "lai", "loss"), "tif$", full.names=TRUE))

lai.con <- file.path(mdl.dir, "lai", "lai_250m_model_data.duckdb")
  con <- dbConnect(duckdb::duckdb(), lai.con, read_only=TRUE)
    db.table <- dbListTables(con)[1]
    dat <- dbGetQuery(con, paste0("SELECT * FROM ", db.table))     
  dbDisconnect(con, shutdown = TRUE)

nidx <- which(dat$NAME %in% "control")
  if(length(nidx) > 0) { dat <- dat[-nidx,] }

#***************************
# Correct loss with no changes
neg.idx <- which(dat$diff >= -0.5 & dat$diff <= 0.5)
  summary(dat$diff[neg.idx])
neg.cells <- dat$cell[neg.idx]

lai.loss <- rast(grep("loss", lai.names, value = TRUE))
  lai.loss[neg.cells] <- 0  

# write raster
loss.out <- file.path(dirname(grep("loss", lai.names, value = TRUE)),
      paste0(rm.ext(basename(grep("loss", lai.names, value = TRUE))), 
	  "_corrected.tif"))
writeRaster(lai.loss, loss,out, overwrite = TRUE) 

plot(lai.loss[[1]], maxcell=5000000, smooth = TRUE, type="interval",  
     breaks = bks, col = hte.mcls)
       title("lai loss effect size")

#***************************
# Correct gain with no changes
pos.idx <- which(dat$diff >= -0.5 & dat$diff <= 0.5)
  summary(dat$diff[pos.idx])
pos.cells <- dat$cell[pos.idx]

lai.gain <- rast(grep("gain", lai.names, value = TRUE))
  lai.gain[pos.cells] <- 0  

# write raster
gain.out <- file.path(dirname(grep("gain", lai.names, value = TRUE)),
      paste0(rm.ext(basename(grep("gain", lai.names, value = TRUE))), 
	  "_corrected.tif"))
writeRaster(lai.gain, gain.out, overwrite = TRUE) 

plot(lai.gain[[1]], maxcell=5000000, smooth = TRUE, type="interval",  
     breaks = bks, col = hte.mcls)
       title("lai gain effect size")








#***************************
# Invert loss
cat("Inverting loss effect size rasters for", country, "\n")

( lai.name <- list.files(file.path(mdl.dir, "lai", "loss"), "tif$", full.names=TRUE) )
  lai.loss <- raster.invert(rast(lai.name))
    writeRaster(lai.loss, lai.name, overwrite = TRUE)
	
( fcov.name <- list.files(file.path(mdl.dir, "fcov", "loss"), "tif$", full.names=TRUE) )
  fcov.loss <- raster.invert(rast(fcov.name))
    writeRaster(fcov.loss, fcov.name, overwrite = TRUE)

remove(lai.name, lai.loss, fcov.name, fcov.loss)
