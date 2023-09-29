library(sf)
library(terra)
library(spatialEco)
library(dplyr)

root = "C:/evans/Colombia/waterfunds"
setwd(file.path(getwd(), "data"))

st_layers("CuencaVerde.gpkg")$name  

#**************************************
# read boundary and format interventions
bdy <- st_read("CuencaVerde.gpkg", "watersheds") |>
    st_union() |>
        st_sf() |> 
          st_cast(to="POLYGON")
    st_geometry(bdy) <- "geometry"


pa <- st_cast(st_read("CuencaVerde.gpkg", layer="protected_areas"), 
              "POLYGON")
  st_geometry(pa) <- "geometry" 
    pa$intervention <- "Protected Area"
	  names(pa)[2] <- "name"
pa <- pa[,c("name", "intervention")]

treatment <- st_read("CuencaVerde.gpkg", "treatments") 
  treatment$intervention <- "Treatment"
  	names(treatment)[2] <- "name"
treatment <- treatment[,c("name", "intervention", "Forest")]

# control <- st_read("CuencaVerde.gpkg", "treatment_controls") 
#   control$intervention <- "Control"
# control <- control[,c("name", "intervention")]
   
#***************************************
# create reference raster
ref <- rast(ext(bdy), resolution = 10, crs=crs(bdy)) 

#***************************************
# distance based variables
lulc <- rast("ESA_lulc10m.tif")
  lulc.att <- levels(lulc)[[1]]
    levels(lulc) <- NULL

if(!file.exists("dist_development.tif")){
  ddev <- rasterDistance(lulc, 50)  
    writeRaster(ddev, "dist_development.tif")
} else {
  ddev <- rast("dist_development.tif")
} 
if(!file.exists("dist_crop.tif")){
  dcrop <- rasterDistance(lulc, 40)
    writeRaster(dcrop, "dist_crop.tif", overwrite=TRUE)
} else {
  dcrop <- rast("dist_crop.tif")
}  
if(!file.exists("dist_stream.tif")){ 
  streams <- st_read("CuencaVerde.gpkg", layer="rivers") 
  streams <- mask(rasterize(streams, lulc, background=0), vect(bdy))
  dstrm <- rasterDistance(streams, 1)
    writeRaster(dstrm, "dist_stream.tif")
} else {
  dstrm  <- rast("dist_stream.tif")
} 
if(!file.exists("dist_water.tif")){
  lakes <- st_read( "CuencaVerde.gpkg", layer="lakes") 
    lakes <- mask(rasterize(lakes, lulc, background=0), vect(bdy))
  wtr <- ifel(lulc == 80, 1, 0)
  lakes <- ifel(lakes + wtr > 0, 1, 0)
  dwater <- rasterDistance(lakes, 1)
    writeRaster(dwater, "dist_water.tif")
} else {
  dwater <- rast("dist_water.tif")
}
 
dist.vars <- c(ddev, dcrop, dstrm, dwater) 
  names(dist.vars) <- c("DistDev", "DistCrop",  "DistStrm", "DistWater")

#*****************************************************
# Read in raster data  
lulc <- rast("ESA_lulc10m.tif")
clim <- rast("BioClim_pca.tif")
 
# read or vectorize forest landcover
if(!file.exists("landcover_vector.gpkg")){   
  forest <- ifel(lulc == "Tree cover", 1, NA)
    forest <- st_as_sf(as.polygons(forest))
      st_write(forest, "landcover_vector.gpkg", "forest")
} else {
  forest <- st_read("landcover_vector.gpkg", "forest")
}
