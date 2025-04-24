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
library(sf)
library(terra)
library(geodata)

root = "C:/evans/Colombia/waterfunds"
setwd(file.path(root, "data"))

#***************************************
# read province for extent
bdy <- st_read("CuencaVerde.gpkg", "watersheds") |>
    st_union() |>
        st_sf() |> 
          st_cast(to="POLYGON")
    st_geometry(bdy) <- "geometry"
ref <- rast(ext(bdy), resolution=1000, crs=crs(bdy))

#***************************************
# download BIOCLIM variables and reproject
clim <- worldclim_country(country = "Colombia", var="bio", 
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
	    names(pca) <- c("pca1","pca2")
writeRaster(pca, "BioClim_pca.tif", overwrite=TRUE)  



#pca[[1]] <- raster.transformation(pca[[1]], "stretch", 0,255)
#pca[[1]][pc$cell] <- as.numeric(pc[,2])

#pca[[2]][pc$cell] <- as.numeric(pc[,3])
#pca[[2]] <- raster.transformation(pca[[2]], "stretch", 0,255)

#pca[[3]][pc$cell] <- as.numeric(pc[,4])
#pca[[3]] <- raster.transformation(pca[[3]], "stretch", 0,255)

#plotRGB(pca, 1, 2, 3, stretch="hist", smooth=TRUE)

