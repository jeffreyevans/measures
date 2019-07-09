########################## START CODE ##########################
#  This will install the required libraries in the R library directory. 
##############################################################
# Set site library to C:/Program Files/R/R-version/library
.Library.site <- file.path(chartr("\\", "/", R.home()), "library")
.libPaths(file.path(chartr("\\", "/", R.home()), "library"))

# Required Edits - Add conditional - if package not installed install
# Refined package list: 

# Install libraries
 # install.packages(c("sp", "spdep", "raster", "rgeos", "spatstat", "rgdal", "maptools", 
  # "spatialEco", "RANN", "RColorBrewer", "spatial", "Matrix", "nlme", "randomForest", 
  # "SDMTools", "rfUtilities", "vegan","yaImpute", "classInt", "cluster", "e1071", "rgl", "rms",
  # "spgrass6", "spgrass7", "RSAGA"), repos="http://cran.us.r-project.org",
  # dependencies=c("Depends", "Imports", "LinkingTo", "Suggests", "Enhances"))
########################## END CODE ##########################

