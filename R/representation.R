suppressMessages(
  lapply(c("sf", "spatialEco", "terra"), 
			require, character.only = TRUE))

cty <- c("brazil", "bhutan", "canada", "colombia", "costa_rica", "peru")
  out.elu <- "C:/evans/PFP/results/representation_elu.csv"
  out.clim <- "C:/evans/PFP/results/representation_bioclimate.csv"

mask.forest = c(TRUE,FALSE)[1]

for(j in cty) {
  cat("Calculating percentages for", j, "\n")
  setwd(file.path("C:/evans/PFP", j, "data"))
  if(j %in% c("brazil", "colombia", "peru")) {
    pa <- st_read(paste0(j, ".gpkg"), "interventions")
      rm.idx <- grep(paste(c("indigenous","Indigenous"), collapse="|"), pa$DESIG)
        if(length(rm.idx > 0)) pa <- pa[-rm.idx,]  
  } else {
    pa <- st_cast(st_read(paste0(j, ".gpkg"), "protected_areas"), "POLYGON")
      iidx <- grep(paste(c("indigenous","Indigenous"), collapse="|"), pa$DESIG)
  	  if(length(iidx) > 0) pa <- pa[-iidx,]
  }
  if(mask.forest) {
    fm <- rast("forest_250m.tif")
      fm[fm<1] <- NA
    elu <- mask(rast("elu.tif"), fm)
	clim <- mask(rast("elu_temperature_moisture.tif"), fm)
  } else {
    elu <- rast("elu.tif")
	clim <- rast("elu_temperature_moisture.tif")
  }
  
  f <- freq(elu, bylayer=FALSE)
    f <- f[-grep("converted", f$value),]
      f <- tapply(f$count, f$value, sum)
        f <- data.frame(ecosystem=names(f), country_count = as.numeric(f))
          f$country_pct <- f$country_count / sum(f$country_count)   
  elu <- mask(elu, pa) 
    f.pa <- freq(elu)
	  f.pa <- f.pa[-grep("converted", f.pa$value),]
        f.pa <- tapply(f.pa$count, f.pa$value, sum)
          f.pa <- data.frame(ecosystem=names(f.pa), pa_count = as.numeric(f.pa))
            f.pa$pa_pct <- f.pa$pa_count / sum(f.pa$pa_count) 
  f <- merge(f, f.pa, by="ecosystem", all=TRUE)	 
    f <- data.frame(country=j, f, fmasked=mask.forest)
	  f[is.na(f)] <- 0   
  write.table(f[order(f$country_count, decreasing = TRUE),], out.elu, sep = ",", 
              row.names = FALSE, append = TRUE, 
              col.names = !file.exists(out.elu)) 

  f <- freq(clim, bylayer=FALSE)
    f <- f[-grep("converted", f$value),]
      f <- tapply(f$count, f$value, sum)
        f <- data.frame(ecosystem=names(f), country_count = as.numeric(f))
          f$country_pct <- f$country_count / sum(f$country_count)   
  clim <- mask(clim, pa) 
    f.pa <- freq(clim)
	  f.pa <- f.pa[-grep("converted", f.pa$value),]
        f.pa <- tapply(f.pa$count, f.pa$value, sum)
          f.pa <- data.frame(bioclimate=names(f.pa), pa_count = as.numeric(f.pa))
            f.pa$pa_pct <- f.pa$pa_count / sum(f.pa$pa_count) 
  f <- merge(f, f.pa, by="bioclimate", all=TRUE)	 
    f <- data.frame(country=j, f, fmasked=mask.forest)
	  f[is.na(f)] <- 0
  write.table(f[order(f$country_count, decreasing = TRUE),], out.clim, sep = ",", 
              row.names = FALSE, append = TRUE, 
              col.names = !file.exists(out.clim)) 
}


#  #**************************************
#  # Ecological representation Function 
#  #   x  raster, classes to be evaluated
#  #   y  polygons (intervention units)
#  #   e  polygon (study area eg., country)
#  #   w  weights (names corresponding with raster classes) 
#  eco.rep(x, y, e, w = NULL) {
#    st_erase <- function(x, y) sf::st_difference(x, sf::st_union(sf::st_combine(y)))
#    d <- st_erase(e, y)
#    gf <- terra::freq(terra::mask(x, terra::vect(d)))
#      gf <- tapply(gf$count, gf$value, sum)
#  	  gf <- gf / sum(gf)	  
#    lf <- terra::freq(terra::mask(x, terra::vect(y)))
#      lf <- tapply(lf$count, lf$value, sum)
#  	  lf <- lf / sum(lf)
#    gf <- data.frame(value=names(gf), pct=as.numeric(gf))
#    lf <- data.frame(value=names(lf), pct=as.numeric(lf))
#    f <- merge(gf, lf, by = "value", all=TRUE)
#      names(f)[2:3] <- c("global.pct", "intervention.pct")
#  	  f[is.na(f)] <- 0
#    return(f)
#  }

