############################################################################
# Code for downloading Data;
#
# 	* MODIS
# 		* Multispectral & QC Flags (Current)
# 		* LAI & FPAR & QC Flags (multi-temporal)
# 		* Landcover & QC Flags (multi-temporal)
# 	* DEM (90m SRTM, static)
# 	* Climate 30 yr normals (PRISM or DAYMET for CONUS, WorldClim for global)
############################################################################
library(raster)
library(spatialEco)
library(sp)
library(ncdf4)
library(rgdal)
library(elevatr)

conus.aea = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
global.aea = "+proj=aea +lat_1=29.5 +lat_2=42.5"
wgs84 = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs "

###############################################
# User defined parameters
###############################################

# Root project directory
root.dir = "C:/..."

# Study area shapefile polygon boundary  
#   (with defined proj), no .shp ext
study.area = "boundary"

# Directory for storing raw data 
raw.dir = "raw"
  dir.create(file.path(root.dir, raw.dir), showWarnings = FALSE)
    raw.dir = file.path(root.dir, raw.dir) 

###############################################
# read data and set projection extent
###############################################

# Read study area polygons
study.area <- readOGR(root.dir, study.area)
proj <- sp::proj4string(study.area)
  if(proj == NA) { warning("Projection must be defined for study area") }

# Create projected extent (e) and WGS84 geographic extent (geo.ext)
e <- raster::extent(study.area)
  e.poly <- as(e, "SpatialPolygons")
    sp::proj4string(e.poly) <- sp::proj4string(study.area)
geo.ext <- extent(sp::spTransform(e.poly, CRS("+init=epsg:4326")))
 
###############################################
# Download elevation DEM/SRTM
############################################### 

## elevatr approach (preferred)
# elevation <- elevatr::get_elev_raster(as(e, "SpatialPolygons"), z = 8) 

##  raster getData approach
#  dir.create(file.path(raw.dir, "elev"), showWarnings = FALSE)
#    elev.dir = file.path(raw.dir, "elev")
#    
#  srtm.tiles <- raster(nrows = 24, ncols = 72, xmn = -180,  
#                       xmx = 180, ymn = -60, ymx = 60)
#    srtm.tiles <- raster::crop(strm.tiles, geo.ext)
#    srtm.tiles <- rasterToPoints(srtm.tiles, spatial=FALSE)[,1:2] 
#      for(i in 1:nrow(srtm.tiles)){
#  	  raster:::.SRTM(lon=srtm.tiles[i,][1], lat=srtm.tiles[i,][2], 
#  	                 download = TRUE, path=elev.dir)   
#  	}
  
###############################################
# Download climate data
###############################################
climate.source = c("prism", "prism.norms", "daymet", "worlclim") 
  if(climate.source == "prism.norms") {
    cat("Downloading 800m Prism 30yr normals", "\n")
    } else if(climate.source == "prism") {
	  message("Please download normals")
    } else if(climate.source == "daymet") {
	  message("Please download normals")
    } else if(climate.source == "worldclim") {
	  cat("Downloading 30 arc-second Worldclim 30yr (1970-2000) normals", "\n")
	  # getData('worldclim', var='tmin', res=0.5, lon=5, lat=45)	  
    } else {
      warning("Not a valid option for climate data")
  }
