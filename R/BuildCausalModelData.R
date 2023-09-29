library(sf)
library(terra)
library(spatialEco)
library(dplyr)

root = "C:/evans/Colombia/waterfunds"
setwd(file.path(root, "data"))

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
st_geometry(pa) <- "geometry"

treatment <- st_read("CuencaVerde.gpkg", "treatments") 
  treatment$intervention <- "Treatment"
  	names(treatment)[2] <- "name"
treatment <- treatment[,c("name", "intervention")]
st_geometry(treatment) <- "geometry"

controls <- st_read("CuencaVerde.gpkg", "treatment_controls") 
  controls$intervention <- "Control"
    controls <- controls[,c("name", "intervention")]
st_geometry(controls) <- "geometry"

#***************************************
# distance based variables
lulc <- rast("ESA_lulc10m.tif")
clim <- rast("BioClim_pca.tif")

dist.vars <- rast(c("dist_development.tif", "dist_crop.tif", 
                    "dist_stream.tif", "dist_water.tif")) 
  names(dist.vars) <- c("DistDev", "DistCrop",  "DistStrm", "DistWater")

#***************************************
# fractional forest
forest <- ifel(lulc == 10, 1, 0)

pct.forest <- extract(forest, pa)
pct.forest <- lapply(unique(pct.forest$ID), function(i) {
   v <- factor(na.omit(pct.forest[which(pct.forest$ID == i),][,2]), 
               levels=c("0","1"))
   prop.table(table(v)) })
pct.forest <- do.call(rbind,pct.forest)
pa$pct.forest <- pct.forest[,2] 

pct.forest <- extract(forest, treatment)
pct.forest <- lapply(unique(pct.forest$ID), function(i) {
   v <- factor(na.omit(pct.forest[which(pct.forest$ID == i),][,2]), 
               levels=c("0","1"))
   prop.table(table(v)) })
pct.forest <- do.call(rbind,pct.forest)
treatment$pct.forest <- pct.forest[,2] 

pct.forest <- extract(forest, controls)
pct.forest <- lapply(unique(pct.forest$ID), function(i) {
   v <- factor(na.omit(pct.forest[which(pct.forest$ID == i),][,2]), 
               levels=c("0","1"))
   prop.table(table(v)) })
pct.forest <- do.call(rbind,pct.forest)
controls$pct.forest <- pct.forest[,2] 

plot(pa["pct.forest"], border=NA)
plot(treatment["pct.forest"], border=NA)
plot(controls["pct.forest"], border=NA)

#***************************************
# combine treatments, controls and pa's 
#interventions <- rbind(treatment, controls)
#  interventions <- rbind(interventions, pa)
#plot(interventions["pct.forest"], border=NA)

#plot(st_geometry(pa), col="blue", border=NA)
#  plot(st_geometry(controls), col="black", border=NA, add=TRUE)
#plot(st_geometry(treatment), col="red", border=NA, add=TRUE)

#***************************************
# Coerce to sf and extract data
mref <- rast("LAI_300m_trend.tif", lyr=1)
  mref[] <- 1
    mref <- mask(mref, vect(bdy))

metric.pts <- as.data.frame(mref, xy=TRUE, cells=TRUE, na.rm=TRUE)
  metric.pts <- st_as_sf(metric.pts, coords = c("x", "y"), 
                         crs = st_crs(bdy), agr = "constant")
    metric.pts <- metric.pts[,-2]

# nrow(metric.pts) == length(unique((metric.pts$cell)))

# Attribute PA's
metric.pts$protected.area <- "no"
i = st_intersects(metric.pts, pa, sparse=TRUE)
  idx <- which(unlist(lapply(i, length))>0)  
metric.pts$protected.area[idx] <- "Protected Area"

# Attribute interventions
metric.pts$intervention <- "no"
i = st_intersects(metric.pts, treatment, sparse=TRUE)
  idx <- which(unlist(lapply(i, length))>0)  
metric.pts$intervention[idx] <- "intervention"

# Attribute controls
metric.pts$control <- "no"
i = st_intersects(metric.pts, controls, sparse=TRUE)
  idx <- which(unlist(lapply(i, length))>0)  
metric.pts$control[idx] <- "control"

# Attribute Forest
forest <- ifel(lulc == 10, 1, 0)
pcov <- function(x) {
  if(length(na.omit(as.numeric(x))) > 1) {
    x <- factor(as.numeric(x), levels=c("0","1"))
    p <- as.numeric(prop.table(table(x)))
  } else {
    p <- rep(NA,2)
  }
  return(p)
}
pforest <- focal(forest, w=5, pcov)
  writeRaster(pforest, "percent_forest.tif")
metric.pts$pct.forest <- extract(pforest, metric.pts)[,2]
		   
# Assign raster values
metric.pts$lulc <- terra::extract(lulc, vect(metric.pts))[,2]
metric.pts$elu <- terra::extract(elu, vect(metric.pts))[,2] 
clim.vals <- terra::extract(clim, vect(metric.pts))
   metric.pts <- cbind(metric.pts, clim.vals[,2:4])
dist.vals <- terra::extract(dist.vars, vect(metric.pts))
  metric.pts <- cbind(metric.pts, dist.vals)
    remove(dist.vals)
	
metric.pts <- na.omit(metric.pts)
  metric.pts <- metric.pts[,-c(8)]  

if(nrow(metric.pts) == length(unique((metric.pts$cell)))) {
  st_write(metric.pts, "model_300m_data.gpkg", delete_layer = TRUE)
} else {
  message("Duplicate cell ID's found")
}
