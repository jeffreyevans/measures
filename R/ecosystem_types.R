# Temp_MoistureClassName
#   Temp_ClassName 
#   Moisture_ClassName
# 
# LF_LC_ClassName
#   LF_ClassName              
#   LC_ClassName
#     Landcover_Type
library(terra)
library(sf)

setwd("C:/evans/PFP/costa_rica/data")
  elu.dir <- file.path(getwd(), "elu")
  dir.create(elu.dir, showWarnings = FALSE)
  global.dir <- "C:/evans/GIS/ECOREGIONS/USGS_ELU"

eco.levels <- c("Realm_ClassName",	"LF_ClassName",	
                "LC_ClassName",	"Temp_ClassName",	
                "Moisture_ClassName", "Temp_MoistureClassName",	
                "World_Ecosystem", "Realm_World_Ecosystem",	
                "Landcover_Type")

if(!exists("elu_30m.tif")) {
  elu <- rast(file.path(global.dir, "USGS_ELU_Global.tif"))
  elu.att <- read.csv(file.path(global.dir, "USGS_ELU.csv"))
  # write subset and crosswalk code
} else {
  elu <- rast("elu_30m.tif")
  elu.att <- read.csv("ELU_attributes.csv")
}

if(inherits(levels(elu)[[1]], "data.frame")){
  v <- levels(elu)[[1]]$Value 
} else {
  v <- sort(unique(elu)[,1])
}

# elu.att <- elu.att[,c("Value",eco.levels)]
elu.att <- elu.att[which(elu.att$Value %in% v),]
elu.classes <- levels(elu) <- elu.att[,c("Value", "World_Ecosystem")]
levels(elu) <- data.frame(Value = elu.att$Value, ELU= paste0(elu.att$Temp_MoistureClassName, 
                          " ", elu.att$LF_ClassName))

writeRaster(elu, "elu_level_02.tif")

Polar Moist Mountains          
Cool Temperate Moist Mountains 
Warm Temperate Moist Mountains 
Sub Tropical Moist Mountains  
Sub Tropical Moist Tablelands  
Tropical Moist Plains          
Tropical Moist Hills           
Tropical Moist Mountains      
Tropical Moist Tablelands      
Sub Tropical Moist Plains  