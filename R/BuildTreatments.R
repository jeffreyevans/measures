# Build treatment/control data
library(sf)
library(spatialEco)

root = "C:/evans/Colombia/waterfunds/data/Cuencas_Ref_Control_Impacto"
  setwd(root)
dat.dir <- "C:/evans/Colombia/waterfunds/data"
  dataset <- "CuencaVerde.gpkg"  
  
type = c("Control", "Impact", "Reference") 

treat <- st_cast(st_read(file.path(dat.dir, "Interv_cuencas_m.shp")), "POLYGON")
  treat <- treat[,c("FID_Interv", "NOMBRE_PRE", "ACCIONES", "Tipo_Resta", "Bosque")] 
    names(treat) <- c("ID", "Name", "Actions", "Type", "Forest", "geometry")
      st_write(treat, file.path(dat.dir, dataset), "treatments")
pa <- st_read(file.path(dat.dir, "AreasProtegidas_CuencaVerde.shp"))
  pa <- st_transform(pa, st_crs(treat))
    st_write(pa, file.path(dat.dir, dataset), "protected_areas")

# Parse, combine and write controls
control <- lapply(list.files(file.path(root, type[1]), "shp$", full.names=TRUE), st_read)
  control <- do.call(rbind, control)
   snames <- rm.ext(list.files(file.path(root, type[1]), "shp$"))
     control$intervention <- type[1]
     control$name <- snames
impact <- lapply(list.files(file.path(root, type[2]), "shp$", full.names=TRUE), st_read)
  impact <- do.call(rbind, impact)
   snames <- rm.ext(list.files(file.path(root, type[2]), "shp$"))
     impact$intervention <- type[2]
     impact$name <- snames
reference <- lapply(list.files(file.path(root, type[3]), "shp$", full.names=TRUE), st_read)
  reference <- do.call(rbind, reference)
   snames <- rm.ext(list.files(file.path(root, type[3]), "shp$"))
     reference$intervention <- type[3]
     reference$name <- snames	 
units <- rbind(control, impact)
  units <- rbind(units, reference)	 
st_write(units, file.path(dat.dir, dataset), "treatment_controls")

# Other data
watersheds <- st_read(file.path(dat.dir, "Cuencas.shp"))
  st_write(watersheds, file.path(dat.dir, dataset), "watersheds")
subwatersheds <- st_read(file.path(dat.dir, "Subcuencas.shp"))
  st_write(subwatersheds, file.path(dat.dir, dataset), "subwatersheds")
municipalities <- st_read(file.path(dat.dir, "Municipios.shp"))
  st_write(municipalities, file.path(dat.dir, dataset), "municipalities")

lulc2018 <- st_read(file.path(dat.dir, "Coberturas_2018.shp"))
  lulc2018 <- lulc2018[,c("nivel_1", "nivel_2", "nivel_3", "cambio", "geometry")]
    names(lulc2018) <- c("class01", "class02", "class03", "change", "geometry") 
lulc2012 <- st_transform(st_read(file.path(dat.dir, 
              "CuencaVerde_Coberturas_2010-2012.shp")), st_crs(lulc2018))
  lulc2012 <- lulc2012 <- lulc2012[,c("LEYENDA3N", "geometry")]
    names(lulc2012) <- c("class01", "geometry")

st_write(lulc2012, file.path(dat.dir, dataset), "LULC2012")
st_write(lulc2018, file.path(dat.dir, dataset), "LULC2018")
