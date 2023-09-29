library(sf)
library(terra)
library(spatialEco)

root = "C:/evans/Colombia/waterfunds"
setwd(root)
dat.dir = file.path(getwd(), "data")

wc.dir <-"C:/evans/ESA/WorldCover/ESA_WorldCover_10m_2021_v200_60deg_macrotile_S30W120" 
imperv.dir <- "C:/evans/ESA/impervious/10m_2021"

r <- rast(file.path(dat.dir, "lai_300m.tif"),lyrs=1)
  bb <- bbox_poly(r)
  bb.geo <- st_transform(bb, st_crs(4326))

#***************************************
# read province for extent
bdy <- st_read(file.path(dat.dir, "CuencaVerde.gpkg"), 
               "watersheds") |>
    st_union() |>
        st_sf() |> 
          st_cast(to="POLYGON")
    st_geometry(bdy) <- "geometry"
ref <- rast(ext(bdy), resolution=10, crs=crs(r))

#*****************************************************
# identify tiles, mosaic and reproject WorldCover
tiles <- st_transform(st_read("C:/evans/ESA/WorldCover/tile_grid.shp"),
                      st_crs(crs(r)))

ids <- st_intersection(bb, tiles)$ID
  f <- list.files(wc.dir, "tif$", full.names=TRUE)
    f <- f[grep(paste(ids,collapse="|"), f)] 

r <- mosaic(sprc(lapply(f,rast)))
  r <- mask(crop(r, ext(bb.geo)), vect(bb.geo))
    lulc <- project(r, ref, method="near")
      lulc <- mask(lulc, vect(bdy))
att <- read.csv("C:/evans/ESA/WorldCover/lulc_classes.csv")	
  levels(lulc) <- att	
	
writeRaster(lulc, file.path(dat.dir, "ESA_lulc10m.tif"))

# # subset lulc classes
# v <- levels(lulc)[,1]
# nv <- v
# nv[c(4,5,6,7)] <- NA
# natural <- subst(lulc, v, nv, raw=TRUE)
# writeRaster(natural, file.path(dat.dir, "natural10m.tif")) 
# 
# nv <- v
# nv[c(1:3,6:10)] <- NA
# dev <- subst(lulc, v, nv, raw=TRUE) 
# writeRaster(dev, file.path(dat.dir, "developed10m.tif")) 


#*****************************************************
# identify tiles, mosaic and reproject impervious
tiles <- st_read("C:/evans/ESA/impervious/10m_2021/tile_grid.shp")

ids <- st_intersection(bb.geo, tiles)$ID
  f <- list.files(imperv.dir, "tif$", full.names=TRUE)
    f <- f[grep(paste(ids,collapse="|"), f)] 

r <- mosaic(sprc(lapply(f,rast)))
  r <- mask(crop(r, ext(bb.geo)), vect(bb.geo))
    r <- project(r, ref, method="near")
writeRaster(r, file.path(dat.dir, "impervios10m.tif") 

#*****************************************************
# Subset and reproject USGS Ecological Land Units (realms) 
ref <- rast(ext(bb), resolution=30, crs=crs(r))
eco <- rast("C:/evans/GIS/ECOREGIONS/USGS_ELU/USGS_ELU_Global.tif")
  eco <- crop(eco, ext(bb.geo)) 
    eco <- project(eco, ref, method="near")
      eco <- mask(eco, vect(bdy))

atts <- read.csv("C:/evans/GIS/ECOREGIONS/USGS_ELU/USGS_ELU.csv")
  ids <- sort(unique(eco)[,1])
    idx <- which(atts$Value %in% ids)
      atts <- atts[idx,]

elu.rc <- paste0(atts$Temp_MoistureClassName, " ", atts$LF_ClassName)
  elu.rc[which(atts$Landcover_Type == "Converted")] <- "converted"
    elu.rc[which(atts$LC_ClassName == "Sparsley or Non vegetated")] <- "non vegetated" 
atts$elu <- elu.rc 

levels(eco) <- atts[,c("Value","elu")]

writeRaster(eco, file.path(getwd(), "data", "elu.tif"), overwrite=TRUE)
  write.csv(atts, file.path(getwd(), "data", "ELU_attributes.csv"))

#*****************************************************
# Process World Settlement Footprint (WSF)
ref <- rast(ext(bb), resolution=10, crs=crs(r))

f <- list.files(file.path(getwd(), "data", "WSF"), "tif$", full.names=TRUE)
wsf <- mosaic(rast(f[1]), rast(f[2]))  
  wsf <- crop(wsf, ext(bb.geo)) 
    wsf <- project(wsf, ref, method="near")
      wsf <- mask(wsf, vect(bdy))
wsf[wsf == 1] <- 255
wsf[wsf == 255] <- 1

writeRaster(wsf, file.path(getwd(), "data", "wsf.tif"), overwrite=TRUE)
