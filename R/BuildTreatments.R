# Performs kNN and builds treatment/control data for causal model
#   
#  Builds treatment/control by pulling raster cell-centers
#  and converting them into points. These points are then
#  intersected with interventions (assigning attributes). 
#  Any point intersecting an intervention unit is considered
#  an "intervention" and outside observations are canidate
#  controls.
# 
#  The kNN parameters are then assigned to the points along with 
#  the ELU classes used for stratification. Once the k=1 controls
#  are subset, everything is writtten to a GeoPackage (gpkg)
# 
#  A forest mask is defined by calculating the proportion of 
#    lai >= tropical(3)
#    lai >= temperate(2.5)
#  through the timeseries then creating a binary mask based on 
#  proportion > 0.10. The lai tropical threshold is defined in 
#  Xe et al., (2023); "Estimation of Leaf Area Index in a Typical 
#  Northern Tropical Secondary Monsoon Rainforest by Different 
#  Indirect Methods"
#
#*************************************************
#*************************************************
# set environment
suppressMessages(
  lapply(c("sf", "spatialEco", "terra", "dplyr", "fst", "tsfeatures", 
           "parallel", "snow", "doSNOW", "stringr", "data.table",
		   "duckdb", "DBI"), 
		   require, character.only = TRUE))
terraOptions(tempdir = "C:/temp")

# Country (PFP)
idx = 4  # indicates country and projection 
  country <- c("bhutan", "costa_rica", "canada", "colombia", "peru", "brazil")[idx]
root = file.path("C:/evans/PFP", country)
setwd(file.path(root, "data"))
  mdl.dir = file.path(root, "model")
  dat.dir = file.path(root, "data")

# Metric and resoultion; 
#   lai 300m 2014-2024, lai 250m 2000 - 2021
#   fcov 300m 2014-2024, fcov 500m 2000 - 2021
metric.type <- c("lai", "fcov")[2]           # metric
mres = c("250m", "300m", "500m")[3]          # Resolution of data
  if(mres == "250m") metric.type = "lai"
  if(mres == "500m") metric.type = "fcov"
  
#*************************************************
# Set intervention date
if(!mres == "300m") {
  switch(country, 
    "brazil" = { intervention.date = as.Date("2014-01-10") },  
    "bhutan" = { intervention.date = as.Date("2017-01-10") },  
    "canada" = { intervention.date = as.Date("2006-01-10") },  
    "colombia" = { intervention.date = as.Date("2020-01-10") },  
    "costa_rica" = { intervention.date = as.Date("2010-01-10") },  
    "peru" = { intervention.date = as.Date("2019-01-10") }) 	
} else {
  intervention.date = as.Date("2015-01-10")
}

# Assign numeric resoultion value to nres
switch(mres, "250m" = { nres = 200 }, 
       "300m" = { nres = 300 }, 
	   "500m" = { nres = 500 })  

#*************************************************
# Parameters used in kNN control selection
update.pa = c(TRUE,FALSE)[2]                              # Update and process protected areas
d = c(25000,50000)[1]                                     # distance subset (25 or 50 KM)   
nsamp = 1                                                 # number of controls per intevention obs 
knn.type = c("mv", "mv.dist", "geo.dist", "rand")[1]      # use multivariate kNN or distance based    
stratify = c(TRUE, FALSE)[1]                              # Apply a stratified control sample
update.pa = c(TRUE, FALSE)[2]                             # Update protected areas, if protected_areas does not exists is run
forest.threshold = c(2, 2.5, 1.5, 2.5, 2.5, 2.5)[idx]     # LAI forest threshold
forest.pct = 0.10                                         # fraction in the timeseries that lai >= forest.threshold 
nt = parallel::detectCores()-4                            # sets number of cores for multiprocessing

# Database connections 
ts.db <- file.path(mdl.dir, metric.type, paste0(metric.type, "_",mres, "_timeseries_data.duckdb"))    # timeseries data
tsfe.db <- file.path(mdl.dir, metric.type, paste0(metric.type, "_",mres, "_tsfe.duckdb"))             # feature extraction timeseries
y.db <- file.path(mdl.dir, metric.type, paste0(metric.type, "_",mres, "_response.duckdb"))            # response variables
ct.db <- file.path(mdl.dir, metric.type, paste0(metric.type, "_",mres, "_control_treatment.duckdb"))  # control/treatment data
dat.db <- file.path(mdl.dir, metric.type, paste0(metric.type, "_",mres, "_model_data.duckdb"))        # final model data

# kNN parameters
knn.vars <- c("pca1", "pca2", "pca3", "DistCrop", "DistDev", 
              "DistWater","DistStrm", "DistRds", "DistTown")

#*************************************************
#*************************************************
# Functions
process.ts <- function(j, start = c(2014, 1), f=36) {
  tnames <- c("entropy", "stability", "lumpiness","hurst", "max_level_shift", "time_level_shift", 
    "max_var_shift", "time_var_shift", "max_kl_shift", "time_kl_shift", "crossing_points", "flat_spots",      
     "trend", "spike", "curvature", "e_acf1", "e_acf10", "x_acf1",          
     "x_acf10", "diff1_acf1", "diff1_acf10", "diff2_acf1", "diff2_acf10", "seas_acf1",       
     "x_pacf5", "diff1x_pacf5", "diff2x_pacf5", "seas_pacf", "ARCH.LM", "arch_acf",        
     "garch_acf", "arch_r2", "garch_r2", "facf")  
  heterogeneity <- function(x, resid.corr = TRUE) {
    x.whitened <- as.numeric(na.contiguous(ar(x)$resid))
	 if(resid.corr) x.whitened = x[1:length(x.whitened)] - x.whitened  
    x.archtest <- arch_stat(x.whitened)
    x.archtest <- arch_stat(as.matrix(x.whitened))
    LBstat <- sum(acf(x.whitened^2, lag.max = 12L, plot = FALSE)$acf[-1L]^2)
    garch.fit <- suppressWarnings(tseries::garch(x.whitened, trace = FALSE))
    garch.fit.std <- residuals(garch.fit)
    x.garch.archtest <- arch_stat(garch.fit.std)
    LBstat2 <- NA
    try(LBstat2 <- sum(acf(na.contiguous(garch.fit.std^2), lag.max = 12L, 
        plot = FALSE)$acf[-1L]^2), silent = TRUE)
    output <- c(arch_acf = LBstat, garch_acf = LBstat2, arch_r2 = unname(x.archtest), 
                garch_r2 = unname(x.garch.archtest))
	  if(resid.corr) output <- output / 100			
    return(output)
  }  
  ts.features <- function(x) {  
    tsf <- c(entropy(x), stability(x), lumpiness(x),
      hurst(x), max_level_shift(x), max_var_shift(x),
      max_kl_shift(x), crossing_points(x),
      flat_spots(x), stl_features(x)[c(3,4,6,7,8)],
      acf_features(x), pacf_features(x), arch_stat(x),
      heterogeneity(x), firstmin_ac(x))
    return(tsf)
  } 
  j <- ts(j, start=start, frequency=f)
  ts.metrics <- tryCatch(
    expr = { 
	  ts.features(j) 
	},
      error = function(cond) { 
	    rep(NA, 34)  
    } )
    names(ts.metrics) <- tnames  
  return(round(ts.metrics, 5))
} 

# kNN function
knn.controls <- function(treat, ctl, knn.vars = NULL, nsamp = 1, d = NULL,
                         sel.type = c("mv", "mv.dist", "geo.dist", "rand")) {
  sel.type = sel.type[1]			
  cat("**** Finding control matches", "\n")
    cat("Treatment n =", nrow(treat), " Control n =", nrow(ctl), "\n") 
  if(sel.type == "mv") {  
    if(is.null(knn.vars))
	  stop("Must provide covariates for multivariate control selection")
    cat("**** Using multivariate kNN control selection", "\n")
    knn.ctl <- RANN::nn2(scale(st_drop_geometry(ctl[,knn.vars])),
                         scale(st_drop_geometry(treat[,knn.vars])), 
		                 k=nsamp)
  } else if(sel.type == "dist.mv") {
    if(is.null(knn.vars))
	  stop("Must provide covariates for multivariate control selection")  
    cat("**** Using distance and multivariate kNN control selection", "\n")
    knn.ctl <- RANN::nn2(scale(cbind(st_drop_geometry(ctl[,knn.vars]),
                         st_coordinates(ctl)[,1:2])),
                         scale(cbind(st_drop_geometry(treat[,knn.vars]),
   					     st_coordinates(treat)[,1:2])), 
		                 k=nsamp)
  }  else if(sel.type == "geo.dist") {	
    cat("**** Using distance-based kNN control selection", "\n")
    knn.ctl <- RANN::nn2(st_coordinates(ctl)[,1:2],
                         st_coordinates(treat)[,1:2], 
		                 k=nsamp)	
  }  else if(sel.type == "rand") {
    cat("**** Using random control selection", "\n") 
    knn.idx <- sample(1:nrow(treat), nrow(treat) * nsamp) 
    return(knn.idx)	
  }
  knn.idx <- knn.ctl$nn.idx
    if(is.numeric(d) & sel.type == "rand") {	
	  message("Multivariate distance not relevant for random sample")
    } else if(is.numeric(d)){
	  distances <- knn.ctl$nn.dists
	    distances <- distances / max(distances)
		  distances <- distances < d
	  for(i in 1:ncol(knn.idx)){
        knn.idx[,i][which(distances[,i] == FALSE)] <- NA 
      }	  
	}
    knn.idx <- as.numeric(knn.idx)
  return(knn.idx)
}

trimmed_quantile <- function(x, trim = 1L, use_unique = TRUE, ...) {
  if (trim < 1) {
    trim <- 1L
    warning("Overruling trim less than 1 with trim = 1L")
  }
  if (floor(trim) != trim) {
    trim <- floor(trim)
    warning("Overruling non-integer trim with floor(trim)")
  }
  this_x <- x
  if (use_unique) {
    this_x <- unique(this_x)
  }
  for(i in seq(from = 1L, to = trim, by = 1)) {
    this_x <- this_x[!(this_x %in% range(this_x))]
  }
  stats::quantile(this_x, ...)
}

expected <- function(x, p = 8, trim = 0) {
  x <- x[!is.na(x)]
  if(trim > 0){
    low <- round(length(x) * trim) + 1    
    high <- length(x) + 1 - low
      if (low <= high) 
        x <- sort.int(x, partial = unique(c(low, high)))[low:high]
  }
  weighted.mean(x, x/p)
}

threshold.pct <- function(x, p = 0.70) {
  x <- x[!is.na(x)]
  length(x[x>=p]) / length(x) 
}

ts.vol <- function(x, p = NULL, m = c("trapezoid", "spline", "step")) {
  x <- as.numeric(x)
  vauc <- function(y, x = NULL, method = m[1]) {
    if(is.null(x)) x = 1:length(y)
    idx <- order(x)
    x <- x[idx]
    y <- y[idx]
    switch(match.arg(arg = method, choices = c("trapezoid", "step", "spline")),
      trapezoid = sum((rowMeans(cbind(y[-length(y)], y[-1]))) * (x[-1] - x[-length(x)])),
      step = sum(y[-length(y)] * (x[-1] - x[-length(x)])),
      spline = stats::integrate(stats::splinefun(x, y, method = "natural"), lower = min(x), upper = max(x))$value
    )
  }
  if(is.null(p)) {
    p <- trimmed_quantile(x[!is.na(x)], trim=2, prob=0.75, na.rm = TRUE)
  }	
    x <- x[x >= p]
  return( as.numeric(vauc(x)) )
}

#*************************************************
# Read metric, calculate max LAI and mask to forest
metric <- rast(paste0(toupper(metric.type), "_", mres, "_trend.tif"))
  if(mres == "250m" | mres == "500m") { 
   dates <- as.Date(names(metric))
  } else {
    d <- unlist(lapply(strsplit(names(metric), "_"), function(x) x[2]))
    dates <- as.Date(unlist(lapply(d, function(j) {
        paste( substr(j,1,4), substr(j,5,6), substr(j,7,8), sep="-" )
      })))
  }

# If required, create mask and forest mask
bdy <- st_read(paste0(country, ".gpkg"), "boundary")

if(!file.exists(paste0("mask_",mres,".tif"))) {
    m <- rast(ext(bdy), resolution = res(metric))
      crs(m) <- crs(metric)
	    m[] <- 1
	      m <- mask(m, bdy)
    writeRaster(m, file.path(dat.dir, paste0("mask_",mres,".tif")), 
	              datatype = "INT2U", overwrite=TRUE)  
} else {
  m <- rast(paste0("mask_",mres,".tif"))
}

if(!file.exists(paste0("forest_",mres,".tif"))) {
  if(metric.type != "lai") {
    lai <- rast(paste0("LAI", "_", mres, "_trend.tif")) 
  }
  forest <- app(lai, fun=\(j) { mean(j >= forest.threshold) } ) 
    forest <- mask(ifel(forest >= forest.pct, 1, 0), m)
      writeRaster(forest, file.path(dat.dir, paste0("forest_",mres,".tif")), 
	                datatype = "INT2U", overwrite = TRUE)  
} else {
  forest <- rast(paste0("forest_",mres,".tif"))
    forest <- ifel(forest == 0, NA, forest)
}

#*************************************************
#*************************************************
# Create 
#   1 - PFP and PA intersections
#   2 - Intervention patches raster
#   3 - tiles index polygons
#  ( lyrs <- st_layers("brazil.gpkg")$name )

bdy <- st_read(paste0(country, ".gpkg"), "boundary" )

# If specified, download and process PA update
if(update.pa | !("protected_areas" %in% st_layers(paste0(country, ".gpkg"))$name) ) {
  library(wdpar)
  e <- st_as_sf(as.polygons(ext(bdy)))
    st_crs(e) <- st_crs(bdy)
      e <- st_transform(e, 4326)
  bpa <- wdpa_fetch(country, wait = TRUE,
    download_dir = rappdirs::user_data_dir("wdpar"))
      bpa <- st_intersection(bpa, e) 
        bpa <- wdpa_clean(bpa)
          bpa <- st_transform(bpa, st_crs(bdy))
            bpa <- st_intersection(st_make_valid(bpa), bdy)  
  bpa <- st_as_sf(as.data.frame(bpa), sf_column_name="geometry")
    tidx <- which(bpa$MARINE == "terrestrial") 
      if(length(tidx) > 0) bpa <- bpa[tidx,]
  bpa <- bpa[,c("NAME", "DESIG_ENG", "IUCN_CAT")]
    names(bpa) <- c("NAME", "DESIG", "IUCN", "geometry") 
  na.idx <- unique(c(which(is.na(bpa$IUCN)), which(is.na(bpa$DESIG)))) 
    if(length(na.idx) > 0) bpa <- bpa[-na.idx,]
  bad.geom <- which(st_geometry_type(bpa) == "GEOMETRYCOLLECTION")
    if(length(bad.geom) > 0) bpa <- bpa[-bad.geom,]
  plot(st_geometry(bpa), col="black", border=NA)
  st_write(bpa, paste0(country, ".gpkg"), "protected_areas", append=FALSE)
}

# Read protected areas
if(country %in% c("brazil", "colombia", "peru")) {
  pa <- st_read(paste0(country, ".gpkg"), "protected_areas")
    pa <- pa[pa$DESIG != "Indigenous Area",]
      # Filter out minimum number of pixels in indvidual PA's
      rm.idx <- vector()
        for(i in unique(pa$NAME)){
          s <- pa[pa$NAME == i,]
            a <- extract(forest, s)
      	    a <- tapply(a[,2], a[,1], FUN=\(i) length(na.omit(i)))
                a <- ifelse(sum(a) < 50, TRUE, FALSE)
          if(a) rm.idx <- append(rm.idx, i)
      }
    if(length(a) > 0) pa <- pa[-which(pa$NAME %in% rm.idx),] 
	
  pfp <- st_cast(st_read(paste0(country, ".gpkg"), "PFP"), "POLYGON")
      st_geometry(pfp) <- "geometry"
    if(!any(st_layers(paste0(country, ".gpkg"))$name %in% "pfp_wo_pa")){	  
      pfpc <- pfp |> st_set_precision(1e5) |> st_combine() |> st_make_valid()
        pac <- pa |> st_set_precision(1e5) |> st_combine() |> st_make_valid()
          rm.pa <- st_cast(st_as_sf(st_difference(st_make_valid(pfpc), st_make_valid(pac))), "POLYGON")
            rm.pa$ID <- 1:nrow(rm.pa) 
          a <- units::drop_units(st_area(rm.pa))
        rm.pa <- rm.pa[which(units::drop_units(st_area(rm.pa)) >= 1389055065),]
     # st_write(rm.pa, paste0(country, ".gpkg"), "pfp_wo_pa", append=FALSE)
	} else if(any(st_layers(paste0(country, ".gpkg"))$name %in% "pfp_wo_pa")) {
      pfp <- st_read(paste0(country, ".gpkg"), "pfp_wo_pa")
        pfp$NAME <- "PFP"
        pfp$DESIG <- "PFP"
        pfp$IUCN <- "not reported"
        pfp <- pfp[,-which(names(pfp) %in% "ID")]	
      pa.pfp <- rmapshaper::ms_simplify(rbind(pa, pfp)) 
        d <- pa.pfp %>%
          st_buffer(500) %>% 
            st_union() %>% 
  	        st_as_sf()
        plot(st_geometry(d))
    #st_write(d, paste0(country, ".gpkg"), "intervention_boundaries", append=FALSE)
    }
  interventions <- st_cast(st_read(paste0(country, ".gpkg"), "interventions"), "POLYGON")
    st_geometry(interventions) <- "geometry"  
} else {
  interventions <- st_cast(st_read(paste0(country, ".gpkg"), "protected_areas"), "POLYGON")
    st_geometry(interventions) <- "geometry"
      iidx <- grep("Indigenous", interventions$DESIG)
	    if(length(iidx) > 0) interventions <- interventions[-iidx,]
}

#*************************************************
#*************************************************
# read required vectors and rasters for control selection

# If required, resample ELU from 30m to functional resoultion and classify to 3-level hierarchy
elu.classes <- read.csv("ELU_attributes.csv")
	
if(!file.exists(paste0("elu_", mres, ".tif"))) {
  elu <- rast(file.path(dat.dir, "elu_30m.tif"))
    levels(elu) <- NULL
      elu <- resample(elu, m, method="near", threads = TRUE)
	    elu.vals <- unique(elu)[,1]
	      elu.att <- elu.classes[,c("Value", "World_Ecosystem")]
		    names(elu.att) <- c("value", "elu")
              elu.idx <- which(elu.att$value %in% elu.vals)	  
                levels(elu) <- elu.att[elu.idx,] 
  writeRaster(elu, file.path(dat.dir, paste0("elu_", mres, ".tif")), 
              datatype = "INT2U", overwrite=TRUE)
}

elu <- rast(paste0("elu_", mres, ".tif"))
  levels(elu) <- NULL
   elu.vals <- unique(elu)[,1]
       elu.idx <- which(elu.classes$Value %in% elu.vals)
         elu.classes <- elu.classes[elu.idx,] 
           elu.classes <- elu.classes[,c("Value", "LF_ClassName", "Temp_MoistureClassName")]
	         elu.classes <- data.frame(value = elu.classes[,1], 
	           elu = paste0(elu.classes$LF_ClassName, " ", elu.classes$Temp_MoistureClassName) ) 
    levels(elu) <- elu.classes 

clim <- rast("BioClim_pca.tif")
dist.vars <- rast("distance_variables.tif") 

#*************************************************
#*************************************************
# Coerce to sf and extract data
# coerce raster cell-centers to points
#
#  intervention type = DESIG
#  intervention unit = NAME
#
mpts <- as.data.frame(forest, xy = TRUE, cells = TRUE, na.rm = TRUE)
  mpts <- st_as_sf(mpts, coords = c("x", "y"), crs = st_crs(bdy), agr = "constant")
    mpts$cell <- extract(m, mpts, cells=TRUE)[,"cell"] 
			
  # Intersect with intervention polygons and assign attributes
  inpa <- st_as_sf(terra::intersect(vect(mpts), vect(interventions)))
    inpa <- inpa[,-2]
    if("YEAR" %in% names(inpa)) inpa <- inpa[-which(names(inpa) %in% "YEAR")]
  outpa <- mpts[which(!mpts$cell %in% inpa$cell),] 
      if("NAME" %in% names(inpa)) outpa$NAME <- "control"
      if("IUCN" %in% names(inpa)) outpa$IUCN <- "control"
      if("DESIG" %in% names(inpa)) outpa$DESIG  <- "control"
	  if("PFP" %in% names(inpa)) outpa$PFP  <- "control" 
        outpa <- outpa[,-2]
  mpts <- rbind(inpa, outpa)
  remove(inpa, outpa)
  	   	   
  # Assign raster values
  mpts$elu <- as.character(terra::extract(elu, vect(mpts))[,2])
  mpts <- cbind(mpts, terra::extract(clim, vect(mpts))[,-1])
  mpts <- cbind(mpts, terra::extract(dist.vars, vect(mpts))[,-1])
    mpts <- na.omit(mpts)

#*************************************************
#*************************************************
# Select treatment controls using kNN 
interventions <- mpts[mpts$NAME != "control",]  
controls <- mpts[mpts$NAME == "control",]   
  if(stratify){
    knn.idx <- lapply(unique(interventions$elu), \(i) {
      treat.sub <- interventions[interventions$elu == i,]
      ctl.sub <- controls[which(controls$elu %in% i),]
	  if(nrow(treat.sub) < 5 | nrow(ctl.sub) < 5) {
	    message(i, " does not have enough observsations")
	    return(NULL)
      } else { 	 
        knn.idx <- knn.controls(treat.sub, ctl.sub, 
                                knn.vars = knn.vars, 
                                nsamp = nsamp, 
      			                sel.type = knn.type)
		      ctl.sub <- ctl.sub[knn.idx,]
            ctl.sub$source_cell <- treat.sub$cell 
          treat.sub$source_cell <- 0	  
        rbind(treat.sub, ctl.sub)
	  }
	})
	names(knn.idx) <- unique(interventions$elu)
      knn.idx <- vctrs::list_drop_empty(knn.idx)
        treat.ctl <- do.call(rbind, knn.idx)
  } else {
    knn.idx <- knn.controls(interventions, controls, 
                            knn.vars = knn.vars, 
                            nsamp = nsamp, d = 0.65, 
  						    sel.type = knn.type)
	    controls <- controls[knn.idx,]
      controls$source_cell <- interventions$cell 	  
	treat.ctl <- rbind(interventions, controls)
  }
  rownames(treat.ctl) <- 1:nrow(treat.ctl)
    treat.ctl.db <- data.frame(st_drop_geometry(treat.ctl), 
                               st_coordinates(treat.ctl)[,1:2])
  con <- dbConnect(duckdb::duckdb(dbdir = ct.db), read_only = FALSE)
    if(!DBI::dbExistsTable(con, "treat_control")) {
      duckdb::dbWriteTable(con, "treat_control", treat.ctl.db)
      dbDisconnect(con, shutdown = TRUE)
    } else {
      duckdb::dbAppendTable(con, "treat_control", treat.ctl.db)
      dbDisconnect(con, shutdown = TRUE)
    }
  remove(interventions, controls, knn.idx, mpts, treat.ctl.db)
gc()

#*************************************************
#*************************************************
# Extract timeseries and write to database
if(!exists("treat.ctl")){
  con <- dbConnect(duckdb::duckdb(dbdir = ct.db), read_only = FALSE)
    db.table <- dbListTables(con)[1]
    treat.ctl <- dbGetQuery(con, paste0("SELECT * FROM ", db.table))
  dbDisconnect(con, shutdown = TRUE)
   ets <- data.frame(cell=treat.ctl[,c("cell", "elu")], 
            extract(metric, treat.ctl[,c("X", "Y")])[,-1])  
 } else {
   ets <- data.frame(cell=st_drop_geometry(treat.ctl)[,c("cell", "elu")], 
                     extract(metric, treat.ctl)[,-1])
}
names(ets)[c(1,2)] <- c("cell", "elu") 
dnames <- paste0("D",gsub("-", "", as.character(dates)))
names(ets)[-c(1,2)] <- dnames
	
cat("  Opening database connection and writing table", "\n")
con <- dbConnect(duckdb::duckdb(dbdir = ts.db), read_only = FALSE)
  if(!DBI::dbExistsTable(con, metric.type)) {
    duckdb::dbWriteTable(con, metric.type, ets, append = FALSE)
    dbDisconnect(con, shutdown = TRUE)
  } else {
    duckdb::dbAppendTable(con, metric.type, ets)
    dbDisconnect(con, shutdown = TRUE)
  }

con <- dbConnect(duckdb::duckdb(dbdir = ts.db), read_only = FALSE)
  db.table <- dbListTables(con)[1]
  uelu <- dbGetQuery(con, paste0("SELECT * FROM ", db.table))[,"elu"]
dbDisconnect(con, shutdown = TRUE)

elu.cts <- table(uelu)
uelu <- unique(uelu)

#*************************************************
#*************************************************
# Derive response variable
# current_pct (cpct), current (current), temporal_volume (tsa) 
#*************************************************
#*************************************************

if(!file.exists(ts.db)) # metric_500m_timeseries_data.duckdb
  warning(basename(ts.db), , " ts.db file does not exists")
if(!file.exists(ct.db)) # metric_500m_control_treatment.duckdb
  warning(basename(ct.db), " ct.db file does not exists")

pre.int.idx <- which(dates < intervention.date)   
post.int.idx <- which(dates >= intervention.date)   
min.max.year <- range(as.numeric(format(dates,"%Y"))) 
pre.idx <- unique(sort(c(grep("2001", dates), grep(min.max.year[1], dates))))
post.idx <- unique(sort(c(grep("2014", dates), grep(min.max.year[2], dates))))
metric.max <- max(global(metric, max, na.rm=TRUE)[,1])

#  con <- dbConnect(duckdb::duckdb(dbdir = y.db), read_only = FALSE)
#    db.table <- dbListTables(con)[1]
#    dbGetQuery(con, paste0("DROP TABLE ", db.table))
#  dbDisconnect(con, shutdown = TRUE)

ct=0
  for(e in uelu) {
    ct=ct+1
    con <- dbConnect(duckdb::duckdb(dbdir = y.db), read_only = FALSE)
      if(DBI::dbExistsTable(con, "y")) {
	    n.elu <- nrow(dbGetQuery(con, paste0("SELECT * FROM y WHERE elu LIKE '%", e, "%'")))
        dbDisconnect(con, shutdown = TRUE)			
      } else {
	    n.elu <- 0
	  }
	  if(n.elu > 0) {
        cat("Already processed", e, ct, "of", length(uelu), "proceding to next ELU", "\n") 	  
	    next 
	  }	
    cat("**** Deriving response variables for", e, ct, "of", length(uelu), "\n") 
      flush.console(); Sys.sleep(0.01)
	con <- dbConnect(duckdb::duckdb(dbdir = ts.db), read_only = FALSE)
	  db.table <- dbListTables(con)[1]
      tdat <- dbGetQuery(con, paste0("SELECT * FROM ", db.table, " WHERE elu LIKE '%", e, "%'"))
   dbDisconnect(con, shutdown = TRUE)
    cells <- tdat$cell
    tdat <- tdat[-c(1,2)]
      mmax = max(tdat, na.rm=TRUE)
      onrow = nrow(tdat)
  	  cat("  Number of rows", nrow(tdat), "\n")	  
  	  cat("  Calculating volume", "\n")	  	  
      tsa <- apply(tdat, MARGIN=1, FUN = ts.vol, p=0.1)
        if(inherits(tsa, "list")) {
          tsa[lengths(tsa) == 0] <- NA
            tsa <- unlist(tsa)
        }
  	  cat("  Calculating baselines", "\n")	
        tsa <- tsa / max(tsa, na.rm=TRUE)
          baseline <- apply(tdat[,pre.idx], MARGIN=1, FUN = expected, p = mmax, trim = 0.10)
            baseline[is.nan(baseline)] <- 0.0001
          current <- apply(tdat[,post.idx], MARGIN=1, FUN = expected, p = mmax, trim = 0.10)
            current[is.nan(current)] <- 0.0001
  	  cat("  Calculating responses", "\n")
        # current.pct <- current / (current + baseline)
	     current.pct <- (current / baseline) / (metric.max - baseline) 
          pchg <- (current - baseline) / baseline
            pchg <- ifelse(pchg > 1, 1, 
                      ifelse(pchg < -1, -1, pchg)) 
              na.idx <- unique(c(which(is.infinite(pchg)), which(is.na(pchg)))) 
                if(length(na.idx) > 0) pchg[na.idx] <- 0
                if(length(na.idx) > 0) pchg[na.idx] <- 0
          responses <- data.frame(elu = e, cell = cells, pre = baseline,    
                                  current = current, diff = current - baseline, 
        						    cpct = current.pct, tsa = tsa, pchg = pchg) 
              remove(pchg, current.pct, baseline, current, tdat)
  	cat("  Opening database connection and merging tables", "\n")  
       con <- dbConnect(duckdb::duckdb(dbdir = y.db), read_only = FALSE)
    	  if(!DBI::dbExistsTable(con, "y")) {
            duckdb::dbWriteTable(con, "y", responses)
  	        dbDisconnect(con, shutdown = TRUE)
    	  } else {
  	        duckdb::dbAppendTable(con, "y", responses)
  	        dbDisconnect(con, shutdown = TRUE)
    	  }
           remove(responses)
  	    gc()
        cat("  Finished", as.character(Sys.time()), "\n")
      cat("\n")   
  } # end elu loop
  
#*************************************************
# Time-series Feature Engineering 
# Use feature engineering to produce covariates
#   representing characteristics of the timeseries
# https://cran.r-project.org/web/packages/tsfeatures/vignettes/tsfeatures.html
#

if(!file.exists(ts.db)) # metric_500m_timeseries_data.duckdb
  warning(basename(ts.db), , " ts.db file does not exists")

con <- dbConnect(duckdb::duckdb(dbdir = ts.db), read_only = FALSE)
  db.table <- dbListTables(con)[1]
  felu <- as.data.frame(table(dbGetQuery(con, paste0("SELECT elu FROM ", db.table))[,1]))
  uelu <- unique( dbGetQuery(con, paste0("SELECT elu FROM ", db.table))[,1] )
dbDisconnect(con, shutdown = TRUE)

ct=0
  for(e in uelu){
    ct=ct+1
    con <- dbConnect(duckdb::duckdb(dbdir = tsfe.db), read_only = FALSE)
      if(DBI::dbExistsTable(con, "x")) {
	    db.table <- dbListTables(con)[1]
	      n.elu <- nrow(dbGetQuery(con, paste0("SELECT * FROM ", db.table, " WHERE elu LIKE '%", e, "%'")))
        dbDisconnect(con, shutdown = TRUE)
      } else {
        n.elu = 0
      }	  
	  if(n.elu > 0) {
        cat("Already processed", e, ct, "of", length(uelu), "proceding to next ELU", "\n") 	  
	    next 
	  }
      cat("***************************************************", "\n")
      cat("**** Time-series Feature Engineering for ELU", e, ct, "of", length(uelu), "\n")
      cat("  Started", as.character(Sys.time()), "\n")
        flush.console(); Sys.sleep(0.01)
      con <- dbConnect(duckdb::duckdb(dbdir = ts.db), read_only = FALSE)
        db.table <- dbListTables(con)[1]
        tdat <- dbGetQuery(con, paste0("SELECT * FROM ", db.table, " WHERE elu LIKE '%", e, "%'"))
      dbDisconnect(con, shutdown = TRUE)	
        cells <- tdat$cell
          tdat <- tdat[-c(1,2)]
    	    if(nrow(tdat) < 50000) {
  		  cat("**** Running base apply function", "n =", nrow(tdat),  "for", ct, "of", length(uelu), "\n")
  		    flush.console()	 
              Sys.sleep(0.01)
            tse <- apply(tdat, MARGIN=1, FUN = function(o) { process.ts(na.omit(o)) })
          } else {
            cat("**** Running multi-threaded apply function", "n =", nrow(tdat), "for", ct, "of", length(uelu), "\n")
  		    flush.console()	 
              Sys.sleep(0.01)		  
            cl <- makeCluster(parallel::detectCores()-2, type='SOCK')
              clusterEvalQ(cl, c(library(tsfeatures)))
                clusterExport(cl, c("process.ts"), envir=environment())
              registerDoSNOW(cl)
                tse <- parApply(cl, tdat, MARGIN=1, FUN = function(o) { process.ts(na.omit(o)) }) 
              stopCluster(cl)
            registerDoSEQ()  		
          }	
    	  tse <- data.frame(elu = e, cell = cells, as.data.frame(t(tse)))
         	cat("  Opening database connection and merging tables", "\n")  
            con <- dbConnect(duckdb::duckdb(dbdir = tsfe.db), read_only = FALSE)
        	  if(!DBI::dbExistsTable(con, "x")) {
                duckdb::dbWriteTable(con, "x", tse)
      	        dbDisconnect(con, shutdown = TRUE)
        	  } else {
      	        duckdb::dbAppendTable(con, "x", tse)
      	        dbDisconnect(con, shutdown = TRUE)
        	  }
          cat("  Finished", as.character(Sys.time()), "\n")
        cat("\n")   
      remove(cells, tse)
    gc()	
  } # end elu loop

#*************************************************
#*************************************************
# Build single database
if(!file.exists(y.db)) # metric_500m_response.duckdb
  warning(basename(y.db), " does not exists")
if(!file.exists(ts.db)) # metric_500m_timeseries_data.duckdb
  warning(basename(ts.db), , " ts.db file does not exists")
if(!file.exists(ct.db)) # metric_500m_control_treatment.duckdb
  warning(basename(ct.db), " ct.db file does not exists")
if(!file.exists(tsfe.db)) # metric_500m_tsfe.duckdb
  warning(basename(tsfe.db), " tsfe.db file does not exists")

con <- dbConnect(duckdb::duckdb(dbdir = ts.db), read_only = FALSE)
  db.table <- dbListTables(con)[1]
  felu <- as.data.frame(table(dbGetQuery(con, paste0("SELECT elu FROM ", db.table))[,1]))
  uelu <- unique( dbGetQuery(con, paste0("SELECT elu FROM ", db.table))[,1] )
dbDisconnect(con, shutdown = TRUE)

ct=0
  for(e in uelu){
    ct=ct+1
    con <- dbConnect(duckdb::duckdb(dbdir = dat.db), read_only = FALSE)
      if(!is.na(dbListTables(con)[1])) {
	    n.elu <- nrow(dbGetQuery(con, paste0("SELECT * FROM data WHERE elu LIKE '%", e, "%'")))
        dbDisconnect(con, shutdown = TRUE)			
      } else {
	    n.elu = 0
	  }
	  if(n.elu > 0) {
        cat("Already processed", e, ct, "of", length(uelu), "proceeding to next ELU", "\n") 	  
	    next 
	  }	

    cat("***************************************************", "\n")
    cat("**** Reading", "treatment/control data for ELU", e, ct, "of", length(uelu), "\n")
      cat("  Started", as.character(Sys.time()), "\n")
        flush.console(); Sys.sleep(0.01)
    cat("  Reading treatment/control data for", e, "\n") 	    
      con <- dbConnect(duckdb::duckdb(dbdir = ct.db), read_only = FALSE)
        tcdat <- dbGetQuery(con, paste0("SELECT * FROM treat_control WHERE elu LIKE '%", e, "%'"))
      dbDisconnect(con, shutdown = TRUE)	
    cat("  Reading and merging dependent variables for", e, "\n")   
      con <- dbConnect(duckdb::duckdb(dbdir = y.db), read_only = FALSE)
        db.table <- dbListTables(con)[1]
		response <- dbGetQuery(con, paste0("SELECT * FROM y WHERE elu LIKE '%", e, "%'"))
      dbDisconnect(con, shutdown = TRUE)	
        if(length(which(tcdat$cell == response$cell))  == nrow(response)) {
          tcdat <- data.frame(tcdat[,-which(names(tcdat) %in% knn.vars)], response[,-c(1:2)],
                              tcdat[,which(names(tcdat) %in% knn.vars)])
        } else {
	      cat("Rows in independent data do not match", "\n")  
	    }							  
        remove(response)
    cat("  Reading and merging independent variables for", e, "\n")   	
      con <- dbConnect(duckdb::duckdb(dbdir = tsfe.db), read_only = FALSE)
        tse <- dbGetQuery(con, paste0("SELECT * FROM x WHERE elu LIKE '%", e, "%'"))
      dbDisconnect(con, shutdown = TRUE) 
        if(length(which(tcdat$cell == tse$cell))  == nrow(tse)) {
          tcdat <- data.frame(tcdat, tse[,-c(1:2)])
        } else {
	      cat("Rows in independent data do not match", "\n")  
	    }
        remove(tse)
      # tcdat <- data.frame(tcdat[,1:2], NAME_IDX = as.numeric(factor(tcdat$NAME)), tcdat[,3:ncol(tcdat)])
    cat("size of dB object in RAM", format(object.size(tcdat), units = "auto"), "\n") 
    cat("Writing dB table for", e, "\n")  	
	con <- dbConnect(duckdb::duckdb(dbdir = dat.db), read_only = FALSE)
      if(!DBI::dbExistsTable(con, "data")) {
        duckdb::dbWriteTable(con, "data", tcdat)
        dbDisconnect(con, shutdown = TRUE)
      } else {
        duckdb::dbAppendTable(con, "data", tcdat)
        dbDisconnect(con, shutdown = TRUE)
      }
	remove(tcdat)
  gc()
}
