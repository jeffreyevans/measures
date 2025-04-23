suppressMessages(
  lapply(c("sf", "spatialEco", "terra",  "dplyr", "duckdb", 
           "DBI", "grf", "ggplot2", "data.table"), 
			require, character.only = TRUE))

#***************************
# Set working directory
path = "C:/evans/PFP"
country <- c("brazil", "bhutan", "canada", "colombia", "costa_rica", "peru")[1]

root = file.path(path, country)
setwd(file.path(root, "data"))
  mdl.dir = file.path(root, "model")
  dat.dir = file.path(root, "data")
  results.dir = file.path("C:/evans/PFP/results", country) 

clrs <- c("darkred","darkorange1","azure3","aquamarine3","darkgreen")
bks <- c( seq(-1,-0.005, 0.10), c(-0.01, 0.01), seq(0.10, 1, 0.1)) 
hte.mcls <- colorRampPalette(clrs)(21)     

#******************************************************
#******************************************************
# Correct LAI effect sizes for zero non gains or losses
#******************************************************
#******************************************************
lai.names <- c(list.files(file.path(mdl.dir, "lai", "gain"), "tif$", full.names=TRUE),
               list.files(file.path(mdl.dir, "lai", "loss"), "tif$", full.names=TRUE))

lai.con <- file.path(mdl.dir, "lai", "lai_250m_model_data.duckdb")
  con <- dbConnect(duckdb::duckdb(), lai.con, read_only=TRUE)
    db.table <- dbListTables(con)
      dat <- dbGetQuery(con, paste0("SELECT * FROM ", db.table[1]))
        if(country == "brazil") {
          dat <- rbind(dat, dbGetQuery(con, paste0("SELECT * FROM ", db.table[2])))
        }  
  dbDisconnect(con, shutdown = TRUE)

nidx <- which(dat$NAME %in% "control")
  if(length(nidx) > 0) { dat <- dat[-nidx,] }

#***************************
# Combine and correct no changes
lai.gain <- rast(grep("gain", lai.names, value = TRUE), lyr=1)  
lai.loss <- rast(grep("loss", lai.names, value = TRUE), lyr=1)  
  lai.loss <- raster.invert(lai.loss)
#writeRaster(lai.loss, gsub("_raw", "", grep("loss", lai.names, value = TRUE)),
#            overwrite = TRUE)

change <- rast(grep("gain", lai.names, value = TRUE), lyr=1)  
  change[] <- NA
    names(change) <- "change"

neg.idx <- which(dat$diff <= 0)
  neg.cells <- dat$cell[neg.idx]
change[neg.cells] <- lai.loss[neg.cells]

pos.idx <- which(dat$diff > 0)
  pos.cells <- dat$cell[pos.idx]
change[pos.cells] <- lai.gain[pos.cells]

zero.idx <- which(dat$diff >= -0.001 & dat$diff <= 0.001)
  zero.cells <- dat$cell[zero.idx]
change[zero.cells] <- 0

plot(change, maxcell=5000000, smooth = TRUE, type="interval",  
     breaks = bks, col = hte.mcls)
       title("LAI change effect size") 

#***************************
# write raster
dir.create(file.path(mdl.dir, "lai", "change"), showWarnings = FALSE)
( es.out <- file.path(mdl.dir, "lai", "change", "change_250m_effect_size.tif") )
  writeRaster(change, es.out, overwrite = TRUE) 

# Add to effect size results raster
ftmp <- paste0(tempfile(), ".tif")
  file.rename(file.path(results.dir, "lai_effect_sizes.tif"), ftmp)
    lai.es <- rast(ftmp)
	if("change" %in% names(lai.es)){
	  lai.es[[which(names(lai.es) %in% "change")]] <- change
	} else {
	  lai.es <- c(lai.es, change)
	}  
    if( any(c("gain","loss") %in% names(lai.es)) ) {
      lai.es <- lai.es[[-which(names(lai.es) %in% c("gain","loss"))]]
    }
plot(lai.es, maxcell=5000000, smooth = TRUE, type="interval",  
     breaks = bks, col = hte.mcls)
       title("LAI change effect size") 
	   
writeRaster(lai.es, file.path(results.dir, "lai_effect_sizes.tif"), 
            overwrite = TRUE) 
  file.remove(ftmp)

remove(dat, change, pos.cells, pos.idx, neg.cells, neg.idx, nidx)
  gc()

#******************************************************
#******************************************************
# Correct fCOV effect sizes for zero non gains or losses
#******************************************************
#******************************************************
fcov.names <- c(list.files(file.path(mdl.dir, "fcov", "gain"), "tif$", full.names=TRUE),
               list.files(file.path(mdl.dir, "fcov", "loss"), "tif$", full.names=TRUE))

fcov.con <- file.path(mdl.dir, "fcov", "fcov_500m_model_data.duckdb")
  con <- dbConnect(duckdb::duckdb(), fcov.con, read_only=TRUE)
    db.table <- dbListTables(con)[1]
    dat <- dbGetQuery(con, paste0("SELECT * FROM ", db.table))     
  dbDisconnect(con, shutdown = TRUE)

nidx <- which(dat$NAME %in% "control")
  if(length(nidx) > 0) { dat <- dat[-nidx,] }

#***************************
# Combine and correct no changes
fcov.gain <- rast(grep("gain", fcov.names, value = TRUE), lyr=1)  
fcov.loss <- rast(grep("loss", fcov.names, value = TRUE), lyr=1)
  fcov.loss <- raster.invert(fcov.loss)  
#writeRaster(fcov.loss, gsub("_raw", "", grep("loss", fcov.names, value = TRUE)),
#            overwrite = TRUE)

change <- rast(grep("gain", fcov.names, value = TRUE), lyr=1)  
  change[] <- NA
    names(change) <- "change"

neg.idx <- which(dat$diff <= 0)
  neg.cells <- dat$cell[neg.idx]
change[neg.cells] <- fcov.loss[neg.cells]

pos.idx <- which(dat$diff > 0)
  pos.cells <- dat$cell[pos.idx]
change[pos.cells] <- fcov.gain[pos.cells]

zero.idx <- which(dat$diff >= -0.001 & dat$diff <= 0.001)
  zero.cells <- dat$cell[zero.idx]
change[zero.cells] <- 0

plot(change, maxcell=5000000, smooth = TRUE, type="interval",  
     breaks = bks, col = hte.mcls)
       title("fCOV change effect size")

#***************************	
# write raster
dir.create(file.path(mdl.dir, "fcov", "change"), showWarnings = FALSE)
( es.out <- file.path(mdl.dir, "fcov", "change", "change_500m_effect_size.tif") )
  writeRaster(change, es.out, overwrite = TRUE) 

# Add to effect size results raster
ftmp <- paste0(tempfile(), ".tif")
  file.rename(file.path(results.dir, "fcov_effect_sizes.tif"), ftmp)
    fcov.es <- rast(ftmp)
	if("change" %in% names(fcov.es)){
	  fcov.es[[which(names(fcov.es) %in% "change")]] <- change
	} else {
	  fcov.es <- c(fcov.es, change)
	}  
    if( any(c("gain","loss") %in% names(fcov.es)) ) {
      fcov.es <- fcov.es[[-which(names(fcov.es) %in% c("gain","loss"))]]
    }	
plot(fcov.es, maxcell=5000000, smooth = TRUE, type="interval",  
     breaks = bks, col = hte.mcls)
       title("fCOV effect sizes") 
	   
writeRaster(fcov.es, file.path(results.dir, "fcov_effect_sizes.tif"), overwrite = TRUE) 
  file.remove(ftmp)

remove(dat, change, pos.cells, pos.idx, neg.cells, neg.idx, nidx)
  gc()
