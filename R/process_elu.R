library(terra)
library(sf)

cty <- c("brazil", "bhutan", "canada", "colombia", "costa_rica", "peru")
  elu.full <- rast("C:/evans/GIS/ECOREGIONS/USGS_ELU/USGS_ELU_Global.tif")
    levels(elu.full) <- NULL
  elu.att.full <- read.csv("C:/evans/GIS/ECOREGIONS/USGS_ELU/USGS_ELU.csv")

for(j in cty) {
  cat("Processing elu's for", j, "\n")
  setwd(file.path("C:/evans/PFP", j, "data"))
  bdy <- st_read(paste0(j, ".gpkg"), "boundary")
  if(!file.exists("forest_250m.tif"))	{
    m <- rast(ext(bdy), resolution = 250)
      crs(m) <- crs(bdy)
  	    m[] <- 1
  		  m <- mask(m, bdy)
    } else {
	  m <- rast("forest_250m.tif")
	    m[] <- 1
		  m <- mask(m,bdy)
	}
    e <- st_as_sf(as.polygons(ext(m)))
      st_crs(e) <- st_crs(m)
        e <- st_transform(e, 4326)
    cat("Cropping and projecting ELU's", "\n")
    elu.geo <- crop(elu.full, e)
      elu <- project(elu.geo, m, method="near")
	    elu <- mask(elu, m)
          v <- unique(elu)[,1]
    elu.att <- elu.att.full[which(elu.att.full$Value %in% v),][,c("Value", "Temp_MoistureClassName", 
	                        "World_Ecosystem", "Realm_World_Ecosystem", "Landcover_Type")]
	  elu.att[which(elu.att$Landcover_Type %in% "Converted"),][,2:4] <- "converted"
	    elu.att <- data.frame(Value = elu.att$Value, bioclime=elu.att$Temp_MoistureClassName, 
	                          ecosystem = elu.att$World_Ecosystem, 
			   				  realm_ecosystem = elu.att$Realm_World_Ecosystem)
		write.csv(elu.att, "ELU_attributes.csv", row.names = FALSE)					  
	  levels(elu) <- elu.att[,c("Value", "realm_ecosystem")]
	    temp.moist <- elu 
	      levels(temp.moist) <- elu.att[,c("Value", "bioclime")]   
	writeRaster(elu, "elu.tif", overwrite=TRUE, datatype="INT2U")
	writeRaster(temp.moist, "elu_temperature_moisture.tif", overwrite=TRUE, datatype="INT2U")
}
