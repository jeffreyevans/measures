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
source("C:/evans/GITS/measures/R/FeatureEngineering.R")
# library(DoubleML)

#*************************************************
root = "C:/evans/Colombia/waterfunds"
dat.dir = file.path(root, "data")
mdl.dir <- file.path(root, "model", "causal")
setwd(root)

# Specify time-series metric
metrics <- c("lai", "fcov")[2]

# Define seasonality (eg., dry 12 and 1:5 and wet 6:11)
season = c("none", "dry", "wet")[2] 

# Specify intervention model (specify which intervention to evaluate)
intervention.model <- c("all", "intervention", "Protected Area")[2]

# Specify dependent variable (y) (model will always retain pre intervention condition)
y = c("metric", "pchg", paste0(c("max.", "diff."), metrics))[1]

p = 4                                      # number of candidate controls, selected using kNN 
nsamp = 2                                  # number of counterfactuals, selected using DTW  
nboot = 1001                               # Number of Bootstrap replicates
cdrop = c(2,4,5,14,15)                     # Irrelevant columns to drop in model
knn.idx = 2:9                              # index of columns to retain
sf = 0.25                                  # sample fraction of cross-validation power test
calculate.ci = c(TRUE, FALSE)[2]           # Calculate Confidence Interval 
sample.wts = c("knn", "wts", "random")[1]  # counterfactual sampling method
n.obs = 5                                  # num of obs for deriving median(max) pre/post metric vlaues

#*************************************************
# Functions

# plot heterogeneous treatment effects
plot_htes <- function(cf_preds, ci = FALSE, z = 1.96) {
  out <- ggplot(
    mapping = aes(
      x = rank(cf_preds$predictions), 
      y = cf_preds$predictions
    )
  ) +
    geom_point() +
    labs(x = "Rank", y = "Estimated Treatment Effect") +
    theme_light()
	if(ci) {
    out <- out +
      geom_errorbar(
        mapping = aes(
          ymin = cf_preds$predictions + z * 
		      sqrt(cf_preds$variance.estimates),
          ymax = cf_preds$predictions - z * 
		      sqrt(cf_preds$variance.estimates)
        )
      )
  }	  
  return(out)
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
	  if(max(global(metric, fun="max", na.rm=TRUE)[,1]) > 9)
	    metric[metric > 9] <- 9
    } else if(metrics == "fcov") {
	  if(min(global(metric, fun="min", na.rm=TRUE)[,1]) < 0)
        metric[metric < 0] <- 0
	  if(max(global(metric, fun="max", na.rm=TRUE)[,1]) > 1)
	    metric[metric > 1] <- 1
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
# Create indexes for specified seasonal subset
cat("using", season, "seasonal subset for timeseries response variable", "\n") 
dry.season <- which(months(dates) %in% unique(months(dates))[c(1:5,12)]) 
wet.season <- which(months(dates) %in% unique(months(dates))[6:11]) 
  if(season == "wet") {
    #metric <- metric[[wet.season]]
	dates.sub <- dates[wet.season]	
	  pre.idx <- which(dates <= intervention.date) 
      post.idx <- which( dates > intervention.date) 
  } else if(season == "dry") {
    #metric <- metric[[dry.season]]
	dates.sub <- dates[dry.season]
	  pre.idx <- which(dates.sub <= intervention.date) 
      post.idx <- which(dates.sub > intervention.date)	
  } else { 
    pre.idx <- which(dates <= intervention.date) 
    post.idx <- which( dates > intervention.date)
  }
ts.idx <- c(pre.idx, post.idx)  

#*************************************************
# Create response variables and write associated 
#  rasters
# metric - median of the n max values in the post period 
# min.metric - minimum of timeseries 
# max.metric - maximum of timeseries 
# pre.metric - median of the n max values in the pre period 
# diff.metric - (pre.metric - metric)
# pchg - (metric - pre.metric) / pre.metric 
# mchg - (max.metric - min.metric)
max.metric = apply(tdat, MARGIN=1, FUN=max, na.rm=TRUE)
min.metric = apply(tdat, MARGIN=1, FUN=min, na.rm=TRUE)

# Take median of 5 maximum observations in pre and post period
full.pre.idx <- which(dates %in% dates.sub[pre.idx]) 
  pre.metric = apply(tdat[,full.pre.idx], MARGIN=1, FUN=function(x) {
    median(nth.values(x, N = 5, smallest = FALSE), na.rm=TRUE)
  }) 
full.post.idx <- which(dates %in% dates.sub[post.idx] & year(dates) == "2023")
  m <- apply(tdat[,full.post.idx], MARGIN=1, FUN=function(x) {
    median(nth.values(x, N = 5, smallest = FALSE), na.rm=TRUE)
  })

responses <- data.frame(cell=dat$cell, metric = m, min.metric = min.metric,  
                  max.metric = max.metric, pre.metric = pre.metric, 
				  diff.metric = pre.metric - m, 
                  pchg = (m - pre.metric) / pre.metric,
				  mchg = max.metric-min.metric)
  responses <- round(responses, 4)
 
  # write dependent variable, pre-intervention metric and percent change raster(s) 
  if(!file.exists(file.path(mdl.dir, paste0("ResponseVarible_", metrics, ".tif")))){
    response.raster <- metric[[1:7]]
      response.raster[] <- NA 
        response.raster[responses$cell] <- responses[,2:8]
          names(response.raster) <- names(responses)[2:8]
    writeRaster(response.raster, file.path(mdl.dir, 
          paste0("ResponseVarible_", metrics, ".tif")),
		  overwrite=TRUE) 
  }
  responses <- responses[,c("cell",y,"pre.metric", "pchg")]
    names(responses) <- c("cell", "y", paste0("pre.", metrics), "pct.chg")

remove(max.metric, min.metric, pre.metric, m)
#*************************************************
# Time-series Feature Engineering on controls
# Use feature engineering to produce covariates
#   representing characteristics of the timeseries
# https://cran.r-project.org/web/packages/tsfeatures/vignettes/tsfeatures.html
if(!file.exists(file.path(mdl.dir, paste0("Control_TSFE_", metrics, ".rds")))){
  cf.idx <- which(dat$intervention == "no")
    ctse <- feature.engineering(tdat[cf.idx,], dat$cell[int.idx]) 
      saveRDS(ctse, file = file.path(mdl.dir, 
	          paste0("Control_TSFE_", metrics, ".rds")))
} else {
  ctse <- readRDS(file.path(mdl.dir, paste0("Control_TSFE_", metrics, ".rds")))
}
names(ctse)[1] <- "cell" 

#**************************************
#**************************************
# Run causal impact model
model.data <- list() 
model.lift <- list() 
model.preds <- list()
model.imp <- list()
model.ci <- list()

#**************************************
# start for loop of interventions

# for(j in intervention.model) { 
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
    itse <- feature.engineering(tdat[int.idx,], dat$cell[int.idx]) 
      saveRDS(itse, file = file.path(mdl.dir, paste0(j, "_Intervention_TSFE_", 
	          metrics, ".rds")))
  } else {
    itse <- readRDS(file.path(mdl.dir, paste0(j, "_Intervention_TSFE_", 
	                metrics, ".rds")))
  }
  names(itse)[1] <- "cell" 
  
  controls <- merge(controls, ctse, by="cell")
    controls <- na.omit(controls)  
  interventions <- merge(interventions, itse, by="cell")
    interventions <- na.omit(interventions)
  rm.idx <- which(names(controls) %in% c("pct.forest", "protected.area", "lulc", "elu", "ID", "intervention", "control"))
    controls <- controls[,-rm.idx] 
    interventions <- interventions[,-rm.idx]
  remove(itse)
  
  #*************************************************
  # Select counterfactual controls using kNN and 
  # pre-intervention DTW matching
  # Thresholds kNN at median multivariate distance
  cat("**** Finding control matches using", sample.wts, "\n")
  
  if(sample.wts == "weights") {
    knn.idx <- lapply(1:nrow(interventions), function(i) {
	  treat <- interventions[i,]
        treat <- cbind(st_drop_geometry(treat[,-1]), 
                   matrix(st_coordinates(treat)[,1:2], ncol=2))
		names(treat)[43:44] <- c("X", "Y")	   
	  ctl <- controls[sample(1:nrow(controls), round(nrow(controls) * 0.10, 0)),]
        ctl <- cbind(st_drop_geometry(ctl[,-1]), st_coordinates(ctl)[,1:2])
	  nidx <- RANN::nn2(data=ctl, query=treat, k=nrow(ctl))
	    wts <- as.numeric(nidx$nn.dis[1,])
          wts <- wts/sum(wts)
	  return(sample(1:nrow(ctl), nsamp, prob = wts)) 
	})
	knn.idx <- unique(unlist(knn.idx))
  } else if(sample.wts == "knn") {	
    nidx <- RANN::nn2(scale(cbind(st_drop_geometry(controls[,-1]),
                      st_coordinates(controls)[,1:2])),
                      scale(cbind(st_drop_geometry(interventions[,-1]),
		   			  st_coordinates(interventions)[,1:2])), 
		              k=nsamp)  
        dist.idx <- list()
      dist.idx[[1]] <- FALSE
    for(i in 1:nsamp) {
      dist.idx[[1]] <- c(dist.idx[[1]], nidx$nn.dist[,i] <= quantile(nidx$nn.dist, p=0.50))
    }	  
    dist.idx <- dist.idx[[1]][-1] 
      knn.idx <- as.numeric(nidx$nn.idx)
        knn.idx <- unique(knn.idx[dist.idx])	
  } else if(sample.wts == "random") {
    knn.idx <- sample(1:nrow(controls), nrow(interventions))
  }  
  control.matches <- controls[knn.idx,]
    cat("**** Selected", length(knn.idx), "control matches", "\n") 
  
  #*************************************************
  # combine model data
  interventions <- st_drop_geometry(interventions)
     interventions <- merge(responses, interventions, by="cell")
       interventions <- data.frame(cell=interventions$cell, treatment=1,
                           	      interventions[,-1])	
   control.matches <- st_drop_geometry(control.matches)
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
  # SMD = ( estimated - control ) / pooled(sdev)
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
  X <- mdata[,c(4:ncol(mdata))] 
  Y <- as.numeric(mdata$y)
  TW <- as.numeric(mdata$treatment) 
    cf <- causal_forest(
      X = X,
      Y = Y,
      W = TW,
	  honesty = TRUE,
      honesty.fraction = 0.5,
      num.trees = nboot, 
	  seed = 42,
      tune.parameters = "all" 
    )
  VI <- variable_importance(cf)
    sel.vars.full <- which(VI / mean(VI) > 0.2) 
      #X <- X[,sel.vars.full]

  if(calculate.ci) {
    parms <- cf$tuning.output$params # cross-validated parameters from full model 	
    sample.idx <- sample(1:nrow(mdata), round(nrow(mdata) * sf, 0))
      train.parms <- parms
        train.parms$X <- model.matrix(~ ., data = X[sample.idx,])
          train.parms$Y <- Y[sample.idx] 
          train.parms$W <- TW[sample.idx] 
    train.forest <- do.call(causal_forest, train.parms)
        VI <- variable_importance(train.forest)
          sel.vars.train <- which(VI / mean(VI) > 0.2) 
    eval.parms <- parms
        eval.parms$X <- model.matrix(~ ., data = X[-sample.idx,])
          eval.parms$Y <- Y[-sample.idx]
          eval.parms$W <- TW[-sample.idx] 
    eval.forest <- do.call(causal_forest, eval.parms)  
      VI <- variable_importance(eval.forest)
          sel.vars.eval <- which(VI / mean(VI) > 0.2) 
    xdat <- data.frame(intercept=1, X)
    
    # Compute a prioritization based on estimated treatment effects.
    # -1, the treatment should reduce the risk of an event occurring.
    priority.cate <- -1 * predict(train.forest, xdat[-sample.idx,])$predictions
    
    # Estimate AUTOC on held out data
    rate <- rank_average_treatment_effect(eval.forest, priority.cate)
      # plot(rate) 
    
    # Compute a prioritization based on baseline risk.
    rf.risk <- regression_forest(X[sample.idx[TW[sample.idx] == 0], ], 
                                 Y[sample.idx[TW[sample.idx] == 0]])
      priority.risk <- predict(rf.risk, X[-sample.idx, ])$predictions
      
    # Test if two RATEs are equal.
    rate.diff <- rank_average_treatment_effect(eval.forest, cbind(priority.cate, priority.risk))
      
    # Construct a 95 % confidence interval
    # A significant result suggests that there are HTEs and that the prioritization rule 
    # is effective at stratifying the sample based on them. Conversely, a non-significant 
    # result suggests that either there are no HTEs or the treatment prioritization rule 
    # does not predict them effectively.
    CI95 <- rate.diff$estimate + data.frame(lower = -1.96 * rate.diff$std.err,
              upper = 1.96 * rate.diff$std.err,
                row.names = rate.diff$target)
      CI95[4,] <- c(round(rate$estimate, 2), round(1.96 * rate$std.err, 2))
        rownames(CI95)[4] <- "RATE AUTOC" 
    model.ci[[j]] <- CI95
  }
  
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

    ## Evaluate distributions using Kolmogorov-Smirnov test
	#y.ctl <- mdata[mdata$treatment == 0,]$y
	#y.treat <- mdata[mdata$treatment == 1,]$y 
    #ks.test(y.treat, ecdf(y.ctl))
  
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
      average_treatment_effect(cf, target.sample = "control")[1], 
	  average_treatment_effect(cf, target.sample="overlap")[1])
      names(lift.df) <- c("intervention", "lift", "median error",  
                          "CATE", "CATT", "CATC", "CATO")
    model.lift[[j]] <- lift.df 

  ## lift is the mean treatment minus mean control estimates
  #cf.means <- tapply(mdata$predictions, mdata$treatment, mean)
  #  cf.lift <- cf.means[2] - cf.means[1] 

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
    rname = file.path(mdl.dir, paste0(j, "_", metrics, ".tif"))
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

  #*************************************************
  # Plots
  pdf(file.path(mdl.dir, paste0(metrics, "_plots_", j, ".pdf")),
      height=8.5, width=11)
  
    # Plot distributions of treatment/control of response variable
    mdata$treatment <- factor(mdata$treatment)
	  levels(mdata$treatment) <- c("control", "treatment")
    ct.den <- ggplot(mdata, aes(x=y, fill=treatment)) + 
      geom_density(alpha=0.5)+
        scale_fill_manual(values=c("red","blue")) +
	      ggtitle("Control/Treatment distribution")
    print(ct.den)
    p1 <- ggplot(mdata, aes_string(x = as.name(VI$variable[1]), y = "predictions")) +
      geom_point(colour = "darkblue", size=0.5) +
        geom_smooth(method = "loess") +
          theme_light()
    p2 <- ggplot(mdata, aes_string(x = as.name(VI$variable[2]), y = "predictions")) +
      geom_point(colour = "darkblue", size=0.5) +
        geom_smooth(method = "loess") +
          theme_light()
    p3 <- ggplot(mdata, aes_string(x = as.name(VI$variable[3]), y = "predictions")) +
      geom_point(colour = "darkblue", size=0.5) +
        geom_smooth(method = "loess") +
          theme_light()
    p4 <- ggplot(mdata, aes_string(x = as.name(VI$variable[4]), y = "predictions")) +
      geom_point(colour = "darkblue", size=0.5) +
        geom_smooth(method = "loess") +
          theme_light()
     lf <- cowplot::plot_grid(p1, p2, p3, p4)   
       print(lf)
	   
    ## Plot heterogeneous treatment effects
    hte <- plot_htes(preds, ci = TRUE)
      print(hte)

   ## Plot estimated effect size with confidence intervals
    spreds <- preds[order(preds$predictions),]
      mm <- range(c(spreds$lower,spreds$upper))
      spreds <- preds[order(preds$predictions),]
      plot(1:nrow(spreds), spreds$predictions, type="l", 
           ylim=mm, xlab="", ylab="estimated effect size", 
    	   main="Estimated effect size with confidence interval")
        lines(1:nrow(spreds), spreds$upper, col="grey")
        lines(1:nrow(spreds), spreds$lower, col="grey")

   ## Plot treatment/control timeseires 
    tidx <- which(dat$cell %in% mdata$cell)
      sub.tdat <- data.frame(cell=dat$cell, tdat)
        sub.tdat <- sub.tdat[tidx,]
          treat.idx <- mdata[which(mdata$treatment == 1),]$cell
    	  control.idx <- mdata[which(mdata$treatment == 0),]$cell
    
    #sub.tdat <- sub.tdat[,c(1,full.post.idx+1)]
    #sub.tdat <- sub.tdat[,c(1,full.pre.idx+1)]
    
    treat.idx <- sample(treat.idx, 100)
	control.idx <- sample(control.idx, 100)
    t1 <- smooth.spline(as.numeric(sub.tdat[sub.tdat$cell == treat.idx[1],][2:ncol(sub.tdat)]))
    plot(dates, t1$y, ylim=c(0,1), type="l", main="TNC Treatment/Control timeseries",
	     ylab="Fractional Cover")
      for(i in treat.idx[-1]){
        t1 <- smooth.spline(as.numeric(sub.tdat[sub.tdat$cell == i,][2:ncol(sub.tdat)]))
        lines(dates, t1$y)
      }
      for(i in control.idx){
        t1 <- smooth.spline(as.numeric(sub.tdat[sub.tdat$cell == i,][2:ncol(sub.tdat)]))
        lines(dates, t1$y, col=rgb(1, 0, 0, max = 1, alpha = 0.4))
      }	 
	legend("bottomright", legend=c("treatment","control"), 
           lty=c(1,1), col=c("black", rgb(1, 0, 0, max = 1, alpha = 0.4)))	
	 
  dev.off() 
  #*************************************************
  # Save model
  save.image(file.path(mdl.dir, paste0(metrics, "_Model_", j, ".RData")))
      gc()
	  	  
#**************************************
# } # end for loop
#**************************************
