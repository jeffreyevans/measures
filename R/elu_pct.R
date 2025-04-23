suppressMessages(
  lapply(c("sf", "spatialEco", "terra"), 
			require, character.only = TRUE))

cty <- c("brazil", "bhutan", "canada", "colombia", "costa_rica", "peru")[-1]

out.elu <- "C:/evans/PFP/results/elu_pct.csv"
out.clim <- "C:/evans/PFP/results/climate_gradient_pct.csv"

for(j in cty) {
  cat("Calculating percentages for", j, "\n")
  setwd(file.path("C:/evans/PFP", j, "data"))
  	lai.db = file.path("C:/evans/PFP", j, "model", "lai", "lai_250m_model_data.duckdb")
    fcov.db = file.path("C:/evans/PFP", j, "model", "fcov", "fcov_500m_model_data.duckdb")

 con <- dbConnect(duckdb::duckdb(), lai.db, read_only=TRUE)
      db.table <- dbListTables(con)
      dat <- dbGetQuery(con, paste0("SELECT * FROM ", db.table))
    dbDisconnect(con, shutdown = TRUE)
	dat <- dat[-which(dat$NAME == "control"),]
	  lai <- c(rast(lai.results[1], lyr=c(1,4)),
	           rast(lai.results[2], lyr=trend.idx),
	           rast(lai.results[3], lyr=trend.idx))
          names(lai) <- c("current", "tsa", "pre", "post")
	
	dat$elu <- extract(lai.elu, dat[,c("X", "Y")])[,2]
	  dat <- data.frame(dat[,c("NAME", "DESIG", "IUCN", "elu")], extract(lai, dat[,c("X", "Y")]))
 	




  
  # Read ELU's and interventions, mask to interventions
  elu <- rast("elu_30m.tif")
  if(j %in% c("brazil", "colombia", "peru")) {
    pa <- st_read(paste0(j, ".gpkg"), "interventions")
      rm.idx <- grep(paste(c("indigenous","Indigenous"), collapse="|"), pa$DESIG)
        if(length(rm.idx > 0)) pa <- pa[-rm.idx,]  
  } else {
    pa <- st_cast(st_read(paste0(j, ".gpkg"), "protected_areas"), "POLYGON")
      iidx <- grep(paste(c("indigenous","Indigenous"), collapse="|"), pa$DESIG)
  	  if(length(iidx) > 0) pa <- pa[-iidx,]
  }
  st_geometry(pa) <- "geometry"
  
  elu <- mask(elu, pa) 
    dat <- levels(elu)[[1]] 
	  names(dat)[2] <- "World_Ecosystem"
    elu.att <- read.csv(file.path("C:/evans/GIS/ECOREGIONS/USGS_ELU", "USGS_ELU.csv"))  
      elu.att <- elu.att[grep(paste(unique(dat$World_Ecosystem), collapse="|"), elu.att$World_Ecosystem),]
        elu.att <- elu.att[!duplicated(elu.att$World_Ecosystem),][,c("World_Ecosystem", "Temp_MoistureClassName")]
            names(elu.att)[2] <- "temp_moisture"
    dat <- merge(dat, elu.att, by="World_Ecosystem")
  
  # Write ELU percentagess
  levels(elu) <- dat[,c("value","World_Ecosystem")]
  f <- freq(elu)
    f <- tapply(f$count, f$value, sum)
      f <- data.frame(country=j, names(f), count = as.numeric(f))
        f$pct <- f$count / sum(f$count)
  	      names(f)[2] <- "Ecosystem"
  write.table(f, out.elu, sep = ",", row.names = FALSE, append = TRUE, 
              col.names = !file.exists(out.elu))
  
  # Write climate gradient percentagess
  levels(elu) <- dat[,c("value","temp_moisture")]
  f <- freq(elu)
    f <- tapply(f$count, f$value, sum)
      f <- data.frame(country=j, names(f), count = as.numeric(f))
        f$pct <- f$count / sum(f$count)
  	      names(f)[2] <- "climate_gradient"
  write.table(f, out.clim, sep = ",", row.names = FALSE, append = TRUE, 
              col.names = !file.exists(out.clim)) 
}
