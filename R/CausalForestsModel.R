#*************************************************
# Causal Ensemble Intervention Impact Model 
#
# Copernicus 300m resolution LAI data (2014-current) 
# Uses kNN and DTW for control matching
# Uses feature engineering for converting time-series 
#   to a design matrix
#
# December through February is dry season (12, 1-2)
# March through May is end of dry season (3-5)
# June through August is rainy season (6-8)
# September through November is end of rainy season (9-11)
# 
#  remove(6-9) to mitigate rainy season effects
#*************************************************
suppressMessages(
  lapply(c("sf", "spatialEco", "terra", "dtw", "dplyr", "lubridate", 
         "grf", "rts", "exactextractr", "tsfeatures", "ggplot2", 
		 "rfUtilities"), require, character.only = TRUE))

#*************************************************
# Set environment and model parameters 
metrics <- c("lai", "fcov")[1]

# dry 12 and 1:5 and wet 6:11
season = c("none", "dry", "wet")[2] 

root = "C:/evans/Colombia/waterfunds"
dat.dir = file.path(root, "data")
mdl.dir <- file.path(root, "model", "causal")
setwd(root)

intervention.model <- c("all", "intervention", "Protected Area")[1]

p = 4            # number of candidate controls, selected using kNN 
nsamp = 2        # number of counterfactuals, selected using DTW  
nboot = 1001     # Number of Bootstrap replicates
make.plots = FALSE
y = c("pchg", "max.lai", "pre.lai", "post.lai", "diff.lai")[1]
cdrop = c(2,4,5,14,15)      # Irrelevant columns to drop in model
knn.idx = 2:9

#*************************************************
# Functions
# x <- as.numeric(intervention[1,][8:ncol(intervention)])
heterogeneity <- function(x) {
    x.whitened <- na.contiguous(ar(x)$resid)
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
    return(output)
}

process.ts <- function(j, start = c(2014, 1), f=36) {
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
      ts.metrics <- ts.features(j)  
    names(ts.metrics)[34] <- "facf"  
  return(ts.metrics)
} 

es.summary <- function(x, p=0.0001) {
  neutral <- length(x[x > -p & x < p])
  s <- summary(x)
  ne.idx <- which(x %in% x[x > -p & x < p])
    if(length(ne.idx) > 0) x <- x[-ne.idx]
  return(c( negative = length(x[x < 0]), neutral = neutral,
            positive=length(x[x > 0]), s) )
}

#**************************************
# read study area, protected areas 
# and treatment/control polygons
#
# st_layers(file.path(dat.dir, "CuencaVerde.gpkg"))
bdy <- st_read(file.path(dat.dir, "CuencaVerde.gpkg"), 
               "watersheds") |>
    st_union() |>
        st_sf() |> 
          st_cast(to="POLYGON")
    st_geometry(bdy) <- "geometry"

#**************************************
# read point data with control and treatment
#   intervention n=831 (%forest > 20% n=293) 
#   control n=1641 (%forest > 20% n=915) 
#   protected areas n=7076 (%forest > 20% n=2462)  
# If lai, subset forest > 20%   
cells <- st_read(file.path(dat.dir, "model_300m_data.gpkg"))
  st_geometry(cells) <- "geometry"
  dat <- na.omit(cells)  
    dat$intervention <- ifelse(dat$protected.area != "no", 
	                      "Protected Area", dat$intervention)  
      if(metrics == "lai") {
        dat <- dat[dat$pct.forest >= 0.20,] 
      }

#*************************************************
# read raster time-series and parse dates
# extract time series to model data, fill NA's in 
# the timeseries and then remove > 6 missing values
metric <- rast(file.path(dat.dir, paste0(toupper(metrics), "_300m_trend", ".tif")))
  d <- unlist(lapply(strsplit(names(metric), "_"), function(x) x[2]))
  dates <- as.Date(unlist(lapply(d, function(j) {
      paste( substr(j,1,4), substr(j,5,6), substr(j,7,8), sep="-" )
    })))
	if(metrics == "lai") {
	  if(min(global(metric, fun="min", na.rm=TRUE)[,1]) < 0)
        metric[metric < 0] <- 0
	  if(max(global(metric, fun="max", na.rm=TRUE)[,1]) < 0)		
	    metric[metric > 9] <- 9
    } else if(metrics == "fcov") {
	  if(min(global(metric, fun="min", na.rm=TRUE)[,1]) < 0)
        metric[metric < 0] <- 0
	  if(max(global(metric, fun="max", na.rm=TRUE)[,1]) > 9)		
	    metric[metric > 9] <- 9
    }	

# pre intervention all of 2014
# post intervention (2015-2023) 
intervention.date = dates[37]

# create timeseries data.frame object
tdat <- terra::extract(metric, vect(dat))[,-1]  
  names(tdat) <- gsub("LAI", toupper(metrics), names(tdat))
    na.idx <- as.data.frame(which(is.na(tdat), arr.ind = TRUE))
	if(nrow(na.idx) > 0) {
      na.cts <- as.data.frame(table(na.idx$row))
	      na.cts <- na.cts[na.cts$Freq > 4,]
          na.idx <- as.numeric(as.character(na.cts[which(na.cts$Freq > 10),]$Var1)) 
        if(length(na.idx) > 0) {
	      dat <- dat[-na.idx,]
	      tdat <- tdat[-na.idx,]
	    }
	}

#*************************************************
# Create response variables (based on dry season)
# If specified in season, subset to season	
cat("using", season, "seasonal subset for timeseries response variable", "\n") 
dry.season <- which(months(dates) %in% unique(months(dates))[c(1:5,12)]) 
wet.season <- which(months(dates) %in% unique(months(dates))[6:11]) 
  if(season == "wet") {
    metric <- metric[[wet.season]]
	  pre.idx <- which(dates[wet.season] <= intervention.date) 
      post.idx <- which( dates[wet.season] > intervention.date) 
  } else if(season == "dry") {
    metric <- metric[[dry.season]]
	  pre.idx <- which(dates[dry.season] <= intervention.date) 
      post.idx <- which( dates[dry.season] > intervention.date)	
  } else { 
    pre.idx <- which(dates <= intervention.date) 
    post.idx <- which( dates > intervention.date)
  }
ts.idx <- c(pre.idx, post.idx)  

max.metric = apply(tdat, MARGIN=1, FUN=max, na.rm=TRUE)
min.metric = apply(tdat, MARGIN=1, FUN=min, na.rm=TRUE)
pre.metric = apply(tdat[,pre.idx], MARGIN=1, FUN=median, na.rm=TRUE)
post.metric = apply(tdat[,post.idx], MARGIN=1, FUN=median, na.rm=TRUE)
responses <- data.frame(cell=dat$cell, max.metric = max.metric, 
                  pre.metric = pre.metric, post.metric = post.metric, 
				  diff.metric = pre.metric - post.metric, 
                  pchg = (post.metric - pre.metric) / pre.metric,
				  mchg = max.metric-min.metric)
  responses <- responses[,c("cell",y)]
    names(responses)[2] <- "y"
remove(max.metric, pre.metric, post.metric)

#*************************************************
# Time-series Feature Engineering on controls
# Use feature engineering to produce covariates
#   representing characteristics of the timeseries
# https://cran.r-project.org/web/packages/tsfeatures/vignettes/tsfeatures.html
if(!file.exists(file.path(mdl.dir, paste0("Control_TSFE_", metrics, ".rds")))){
  cf.idx <- which(dat$intervention == "no")
  ctse <- setNames(vector("list", length(cf.idx)), dat$cell[cf.idx])
    for(i in 1:length(cf.idx)) {
	  v <- cf.idx[i]
      tryCatch({
        fe <- process.ts(na.omit(as.numeric(tdat[v,])))
          }, error=function(e){ fe <- rep(NA,34) })
      ctse[[v]] <- fe
    }	
  ctse <- do.call(rbind, ctse)
    ctse <- data.frame(cell=dat$cell[cf.idx], ctse)  
      saveRDS(ctse, file = file.path(mdl.dir, 
	          paste0("Control_TSFE_", metrics, ".rds")))
} else {
  ctse <- readRDS(file.path(mdl.dir, paste0("Control_TSFE_", metrics, ".rds")))
}

#**************************************
#**************************************
# Run causal impact model
model.data <- list() 
model.lift <- list() 
model.preds <- list()
model.imp <- list()

j = intervention.model
 cat("**** Running model for:", j, "**** ", "\n")

  #*************************************************
  # Create control/experimental units
  # Interventions (cat_manejo)
  if(intervention.model == "all") {  
    int.idx <- which(dat$intervention != "no")
      interventions <- dat[int.idx,]
    cf.idx <- which(dat$intervention == "no")
      controls <- dat[cf.idx,]
  } else {
    int.idx <- which(dat$intervention == j)
      interventions <- dat[int.idx,]
    cf.idx <- which(dat$intervention == "no")
      controls <- dat[cf.idx,]
  }
  
  #*************************************************
  # Time-series Feature Engineering on interventions    
  if(!file.exists(file.path(mdl.dir, paste0(j, "_Intervention_TSFE_", metrics, ".rds")))){
  itse <- setNames(vector("list", length(int.idx)), dat$cell[int.idx])
    for(i in 1:length(int.idx)) {
	  v <- int.idx[i]
      tryCatch({
        fe <- process.ts(na.omit(as.numeric(tdat[v,])))
          }, error=function(e){ fe <- rep(NA,34) })
      itse[[v]] <- fe	
    }	 	
    itse <- do.call(rbind, itse)
      itse <- data.frame(cell=dat$cell[int.idx], itse)
        saveRDS(itse, file = file.path(mdl.dir, paste0(j, "_Intervention_TSFE_", 
	            metrics, ".rds")))
  } else {
    itse <- readRDS(file.path(mdl.dir, paste0(j, "_Intervention_TSFE_", 
	                metrics, ".rds")))
  }
  controls <- merge(controls, ctse, by="cell")
    controls <- na.omit(controls)  
  interventions <- merge(interventions, itse, by="cell")
    interventions <- na.omit(interventions)
  rm.idx <- which(names(controls) %in% c("protected.area", "lulc", "elu", "ID", "intervention", "control"))
    controls <- controls[,-rm.idx] 
    interventions <- interventions[,-rm.idx]
  remove(itse)
  #*************************************************
  # Select counterfactual controls using kNN and 
  # pre-intervention DTW matching
  cat("**** Finding control matches", "\n")
  nidx <- RANN::nn2(scale(cbind(st_drop_geometry(controls[,-1]),
                    st_coordinates(controls)[,1:2])),
                    scale(cbind(st_drop_geometry(interventions[,-1]),
					st_coordinates(interventions)[,1:2])), 
		            k=nsamp)$nn.idx
      # nidx <- nidx[,-1]
    control.matches <- controls[unique(as.numeric(nidx)),]
  
  rm(nidx)
  #*************************************************
  # combine model data
  interventions <- st_drop_geometry(interventions[,-knn.idx])
    interventions <- merge(responses, interventions, by="cell")
      interventions <- data.frame(cell=interventions$cell, treatment=1,
                          	      interventions[,-1])	
  control.matches <- st_drop_geometry(controls[,-knn.idx])
    control.matches <- merge(responses, control.matches, by="cell")  
      control.matches <- data.frame(cell=control.matches$cell, treatment=0,
                          	        control.matches[,-1])	
  mdata <- rbind(interventions, control.matches)
 	
  remove(interventions, control.matches)
  #*************************************************
  # Screen for collinear and multi-collinear variables
   all.vars <- names(mdata)[-c(1:3)]    
    cl.vars <- try(collinear(mdata[-c(1:3)], p = 0.85, nonlinear = FALSE, 
	               p.value = 0.001))
	  if(exists("cl.vars")){
        if(length(cl.vars) > 0)
          mdata <- mdata[,-which(names(mdata) %in% cl.vars)]
	  }	
    mc <- multi.collinear(mdata[-c(1:3)], p = 0.05, perm = TRUE, n = 99)
      mc.vars <- mc[which( mc$frequency > 5 ),]$variables 
        if(length(mc.vars) > 0) 
          mdata <- mdata[,-which(names(mdata) %in% mc.vars)]
    cat("The following variables were dropped:",
      all.vars[which(is.na(match(all.vars, names(mdata))))], "\n")

  #*************************************************
  #*************************************************
  # Causal Impact Analysis 
  #
  # Estimate heterogeneous treatment effects as the
  # difference between the expected outcome in the 
  # treatment minus the expected outcome in the control
  # SMD = ( estimated - control ) / pooled(SD)
  #
  # With continuous data the lift is the mean in the treatment minus 
  # the mean in the control so, if the treatment yielded a 5.0 
  # average LAI and the control yielded 3.5, then the lift would 
  # be [5.0 - 3.50 = 1.5]
  #
  # For dichotomous (e.g., treated, non-treated), the lift is the 
  # probability of the desired outcome in the treatment minus the 
  # probability of the desired outcome in the control. 
  # For example, if 55% of the treatment showed increased LAI, 
  # while 50% in the control showed increased LAI, 
  # our lift would be [0.55 - 0.50 = 0.05] 
  # 
  #*************************************************
  # Causal Ensemble model 
  # Heterogeneous treatment effect estimation
  # Adversarial method using Fast Gradient Sign Method (FGSM)
  
  cat("**** Causal Impact Ensemble Model for", j, "\n")
  # effect size regression 
  cf <- causal_forest(
    X = model.matrix(~ ., data = mdata[,c(4:ncol(mdata))]),
    Y = mdata$y,
    W = as.numeric(mdata$treatment),
    num.trees = nboot,
    seed = 42) 
  
  #*************************************************
  # Write model results:  
  cat("**** Writing Model Results for:", j, "\n")
 
  #************************************************* 
  # generate predictions
  preds <- predict( object = cf, estimate.variance = TRUE)
    preds$lower <- preds$predictions - 1.96 * sqrt(preds$variance.estimates)
    preds$upper <- preds$predictions + 1.96 * sqrt(preds$variance.estimates)
    preds$error <- preds$debiased.error + preds$excess.error
  model.preds[[j]] <- preds
  
  mdata <- data.frame(predictions = preds$predictions, mdata)  

  #*************************************************  
  # lift - using upper 50% estimate percentile
  # median.error - using sum of: 
  #   debiased.error - estimates of the R-loss (Nie and Wager 2017). 
  #   excess.error - jackknife estimates of the Monte-carlo error (Wager et al., 2014) 
  #                  representing how unstable estimates are.
  # CATE - conditional average treatment effect on the full sample
  # CATT - conditional average treatment effect on the treated sample
  # CATC - conditional average treatment effect on the control sample
  lift <- lm(y ~ treatment, mdata[preds$predictions > 
             median(preds$predictions),])	  
	lift.df <- data.frame(j, coef(lift)[2], median(preds$error), 
      average_treatment_effect(cf, target.sample = "all")[1],
  	  average_treatment_effect(cf, target.sample = "treated")[1],
      average_treatment_effect(cf, target.sample = "control")[1])  
      names(lift.df) <- c("intervention", "lift", "median error",  
                          "CATE", "CATT", "CATC")
    model.lift[[j]] <- lift.df 

  #*************************************************  
  # Variable importance
  VI <- cf %>% 
    variable_importance() %>% 
      as.data.frame() %>% 
        mutate(variable = colnames(cf$X.orig)) %>% 
          arrange(desc(V1))
  VI <- VI[-nrow(VI),]
    VI[,1] <- VI[,1] / max(VI[,1])
	names(VI)[1] <- "importance"
  model.imp[[j]] <- VI 

  #*************************************************
  # Create spatial predictions (rasters)
  if(exists("preds")){
    pred.idx <- which(mdata$treatment == 1)  
    rname = file.path(mdl.dir, "causal", paste0(j, "_", metrics, ".tif"))
      es.raster <- metric[[1]]
	    es.raster[] <- NA
          names(es.raster) <- "effect.size" 
       es.raster[mdata[pred.idx,]$cell] <- round(preds[pred.idx,]$predictions,6)	
    error.name = file.path(mdl.dir, "causal", 
                paste0(j, "_", metrics, "_errors.tif"))
      error.raster <- metric[[1]]
	    error.raster[] <- NA
          names(error.raster) <- "error"  
      error.raster[mdata[pred.idx,]$cell] <- round(preds[pred.idx,]$error,6)
    writeRaster(c(es.raster, error.raster), rname, datatype="FLT4S", 
	            overwrite=TRUE)  
  }  
  save.image(file.path(mdl.dir, paste0(metrics, "_Model_", j, ".RData")))
      gc()
	  	  
#**************************************
#**************************************
