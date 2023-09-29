suppressMessages(
  lapply(c("sf", "spatialEco", "terra", "sf", "dplyr", "osmdata"), 
		 require, character.only = TRUE))

root = "C:/evans/Colombia/waterfunds"
setwd(root)
dat.dir = file.path(getwd(), "data")

geo.prj = st_crs(4326)

#***************************************
# read protected area polygons, dissolve
# internal borders 

#st_layers(file.path(dat.dir, "CuencaVerde.gpkg"))
bdy <- st_read(file.path(dat.dir, "CuencaVerde.gpkg"), 
               "watersheds") |>
    st_union() |>
        st_sf() |> 
          st_cast(to="POLYGON")
    st_geometry(bdy) <- "geometry"

cv.prj <- st_crs(bdy) 

geo.ext <- st_transform(bbox_poly(bdy), geo.prj)
geo.bb  <- terra::ext(geo.ext)[c(1,3,2,4)]

#***************************************
# Subest HydroSheds data

lakes <- st_make_valid(st_read("C:/evans/GIS/HydroSheds/HydroLAKES.gdb", 
                       "HydroLAKES_polys_v10"))
  lakes <- st_intersection(lakes, geo.ext)
    lakes <- st_intersection(st_transform(lakes, st_crs(bdy)), bdy)  
      st_write(lakes, file.path(dat.dir, "CuencaVerde.gpkg"), "lakes")

rivers <- st_read("C:/evans/GIS/HydroSheds/HydroRIVERS.gdb", 
                  "HydroRIVERS_v10")
  rivers <- st_intersection(rivers, geo.ext)
    rivers <- st_intersection(st_transform(rivers, st_crs(bdy)), bdy)
      st_write(rivers, file.path(dat.dir, "CuencaVerde.gpkg"), "rivers")

#**************************************************
# Query OSM villages, towns and cities
# available_features()
# available_tags("place")
qkey = "place"  
n <- c("osm_id", "name", "place")
villages <- opq(bbox=geo.bb) %>%
    add_osm_feature(key = qkey, value = "village") %>%
	  osmdata_sf()
        villages$osm_points$type <- "village"
        villages <- st_intersection(st_transform(villages$osm_points, 
		                            st_crs(bdy)), bdy)
        villages <- villages[,n]  
		
town <- opq(bbox=geo.bb) %>%
    add_osm_feature(key = qkey, value = "town") %>%
      osmdata_sf()
        town$osm_points$type <- "town"
        town <- st_intersection(st_transform(town$osm_points, 
		                            st_crs(bdy)), bdy)
        town <- town[,n]      
	   
settlements <- rbind(villages, town)
  settlements <- st_difference(settlements)
st_write(settlements, file.path(dat.dir, "CuencaVerde.gpkg"), 
         "towns")

#**************************************************
# Query OSM roads
n <- c("osm_id", "type")
qkey = "highway"  
rds <- opq(bbox=geo.bb) %>%
	     add_osm_feature(key = qkey, value = "primary") %>%
	       osmdata_sf()				
             rds$osm_lines$type <- "primary"
             rds <- st_intersection(st_transform(rds$osm_lines, 
		                            st_crs(bdy)), bdy)
             rds <- rds[,n]  
st_write(rds, file.path(dat.dir, "CuencaVerde.gpkg"), "roads")
