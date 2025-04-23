# Create covariates for kNN control selection 
# and causal model
#
# BioClim_pca - 3 laoding of BioClim PCA 
# dist_crop - distance to crop/ag       
# dist_development - distance to development
# dist_road - distance to roads     
# dist_stream - distance to streams   
# dist_town - distance to towns/settletments        
# dist_water - distance to water bodies 
#
# BioClim (not used, for deriving PCA)        
# elu (used in control stratifiaction)              
# forest_300m, forest_250m (mask for forested areas)      
# lulc10m (not used, for providing landcover classes)           
#
library(sf)
library(terra)
library(spatialEco)
library(dplyr)
library(geodata)

mres = c("250m", "300m")[1]
if(mres == "250m") { out.res = 250 } else { out.res = 300 }

idx = 1  # indicates country 
  country <- c("brazil/ARPA", "bhutan", "canada/great_bear", "costa_rica",  
               "colombia/HECO", "peru")[idx]
    root = file.path("C:/evans/PFP", country) 
      if(length(grep("/", country)) > 0) { 
        pfp = dirname(country) #(country) 
        country = dirname(country) 
      } else { 
        pfp = country 
      }
setwd(file.path(root, "data"))

gp <- paste0(pfp, ".gpkg")

#**************************************
# read reference polygons and raster
# ( lyrs <- st_layers(paste0(pfp, ".gpkg"))$name )
bdy <- st_read(gp, "boundary")

if(!file.exists(paste0("mask", mres, ".tif"))) {
  ref <- rast(vect(bdy), resolution = out.res, crs = crs(bdy))
    ref <- rasterize(bdy, ref)
	  writeRaster(ref, paste0("mask", mres, ".tif"), 
	              datatype="INT1U")
} else {
  ref <- rast(paste0("mask", mres, ".tif"))
}

#***************************************
# Create forest mask (250m or 300m)
if(!file.exists(paste0("forest_", mres, ".tif"))) {
  f <- rast("forest_10m.tif")
    forest <- resample(f, ref, "near") 
      writeRaster(forest, paste0("forest_", mres ,".tif"),
                  datatype="INT1U")
} else {
  forest <- rast(paste0("forest_", mres, ".tif"))
}

#***************************************
# distance based variables
dref <- rast(ext(bdy), resolution = 90, crs=crs(bdy)) 
lulc <- rast("lulc10m.tif")
  levels(lulc) <- NULL
    lulc <- resample(lulc, dref, method = "near")
lulc.vals = unique(lulc)[,1]

if( 40 %in% lulc.vals ) {
  if(!file.exists("dist_crop.tif")){
    dcrop <- rasterDistance(lulc, 40)
      writeRaster(dcrop, "dist_crop.tif", overwrite=TRUE)
  } else {
    dcrop <- rast("dist_crop.tif")
  }  
}
if( 50 %in% lulc.vals ) {
  if(!file.exists("dist_development.tif")){
    ddev <- rasterDistance(lulc, 50)  
      writeRaster(ddev, "dist_development.tif", overwrite=TRUE)
  } else {
    ddev <- rast("dist_development.tif")
  }
}

if(any(st_layers(gp)$name == "lakes")) { 
  if(!file.exists("dist_water.tif")){
    water <- st_read(gp, layer="lakes") 
    water <- mask(rasterize(water, dref, background=0), vect(bdy))  
    dwater <- rasterDistance(water, 1)
      writeRaster(dwater, "dist_water.tif", overwrite=TRUE)
  } else {
    dwater <- rast("dist_water.tif")
  }
}

if(any(st_layers(gp)$name == "streams")) {  
  if(!file.exists("dist_stream.tif")){
    streams <- st_read(gp, layer="streams") 
    streams <- mask(rasterize(streams, dref, background=0), vect(bdy))  
    dstrm <- rasterDistance(streams, 1)
      writeRaster(dstrm, "dist_stream.tif", overwrite=TRUE)
  } else {
    dstrm  <- rast("dist_stream.tif")
  }
}
  
if(any(st_layers(gp)$name == "roads")) {   
  if(!file.exists("dist_road.tif")){ 
    rds <- st_cast(st_read(gp, layer="roads"), "LINESTRING")
    rds <- mask(rasterize(rds, dref, background=0), vect(bdy)) 	
    drds <- rasterDistance(rds, 1)
      writeRaster(drds, "dist_road.tif", overwrite=TRUE)
  } else {
    drds  <- rast("dist_road.tif")
  }
}

if(any(st_layers(gp)$name == "towns")) { 
  if(!file.exists("dist_town.tif")){
    towns <- st_buffer(st_read(gp, layer="towns"), 1000)
    towns <- mask(rasterize(towns, dref, background=0), vect(bdy))  	
    dtown <- rasterDistance(towns, 1)
      writeRaster(dtown, "dist_town.tif", overwrite=TRUE)
  } else {
    dtown  <- rast("dist_town.tif")
  }
}

dist.vars <- c(dcrop, ddev, dwater, dstrm, drds, dtown)
  names(dist.vars) <- c("DistCrop", "DistDev", "DistWater", "DistStrm", "DistRds", "DistTown")
writeRaster(dist.vars, "distance_variables.tif", overwrite=TRUE)

################################################################
# process WorldClim climate data
# Bioclim variable descriptions
# BIO1 = Annual Mean Temperature
# BIO2 = Mean Diurnal Range (Mean of monthly (max temp - min temp))
# BIO3 = Isothermality (BIO2/BIO7) (×100) 
# BIO4 = Temperature Seasonality (standard deviation ×100) 
# BIO5 = Max Temperature of Warmest Month
# BIO6 = Min Temperature of Coldest Month 
# BIO7 = Temperature Annual Range (BIO5-BIO6)
# BIO8 = Mean Temperature of Wettest Quarter
# BIO9 = Mean Temperature of Driest Quarter
# BIO10 = Mean Temperature of Warmest Quarter
# BIO11 = Mean Temperature of Coldest Quarter
# BIO12 = Annual Precipitation
# BIO13 = Precipitation of Wettest Month
# BIO14 = Precipitation of Driest Month
# BIO15 = Precipitation Seasonality (Coefficient of Variation)
# BIO16 = Precipitation of Wettest Quarter
# BIO17 = Precipitation of Driest Quarter
# BIO18 = Precipitation of Warmest Quarter
# BIO19 = Precipitation of Coldest Quarter

ref <- rast(ext(bdy), resolution = 1000, crs=crs(bdy))

#***************************************
# download BIOCLIM variables and reproject
dir.create(file.path(getwd(), "tmp"), 
           showWarnings = FALSE)

clim <- worldclim_country(country = country, var="bio", 
                  path=file.path(getwd(),"tmp"))
  names(clim) <- c("Temp", "DiurnalRange", "Isothermality",
       "TempSeasonality", "MaxTempWarmMt", "MinTempColdMt",
       "TempRange", "TempWetQt", "TempDryQt", "TempWarmQt",
       "TempColdQt", "Precp", "PrecpWetMt", "PrecpDruMt",
       "PrecpSeasonality", "PrecpWetQt", "PrecpDryQt",
       "PrecpWarmQt", "PrecpColdQt")
 
clim <- project(clim, ref, method="bilinear")  
  clim <- mask(clim, vect(bdy)) 
writeRaster(clim, "BioClim.tif", overwrite=TRUE)

#***************************************
# Perform PCA 
rd <- as.data.frame(clim, cells=TRUE, na.rm=TRUE)
pc <- prcomp(rd[,-1], center = TRUE, scale = TRUE)
  summary(pc)
pc <- data.frame(cell=rd$cell, as.data.frame(pc$x[,1:3]))
  pca <- clim[[1:3]]
    pca[[1]][pc$cell] <- as.numeric(pc[,2])
	  pca[[2]][pc$cell] <- as.numeric(pc[,3])
	    names(pca) <- c("pca1","pca2","pca3")
writeRaster(pca, "BioClim_pca.tif", overwrite=TRUE)  

plotRGB(pca, 1, 2, 3, stretch="hist", smooth=TRUE)
