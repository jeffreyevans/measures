##################################################
##################################################
#      Adversarial Causal Ensemble Model 
#
# Copernicus 1000m resolution LAI data 
# Uses kNN and DTW for control matching
# Uses feature engineering for time-series parameters
##################################################
##################################################
# [1] "Protected Areas",
# Berau Community Initiatives:
#   [2] "IND002.3_Social_Forestry_Areas",
#   [3] "IND002.4_Barat_KPH_FMU_Model_Area",
#   [4] "IND002.5_Wehea_Kelay_Ecosystem_Essential_Area"	

# add require packages
suppressMessages(
  lapply(c("sp", "sf", "raster", "spatialEco", "dtw", "dplyr", "lubridate", 
           "grf", "rts", "exactextractr", "tsfeatures", "ggplot2"), 
		   require, character.only = TRUE))

setwd("C:/evans/measures/Indonesia")
  ddir = file.path(getwd(), "data")

################################################
# Functions
# x <- as.numeric(intervention@data[1,][8:ncol(intervention)])
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

process.ts <- function(j) {
  ts.features <- function(x) {
    tsf <- c(entropy(x), stability(x), lumpiness(x),
      hurst(x), max_level_shift(x), max_var_shift(x),
      max_kl_shift(x), crossing_points(x),
      flat_spots(x), stl_features(x)[c(3,4,6,7,8)],
      acf_features(x), pacf_features(x), arch_stat(x),
      heterogeneity(x), firstmin_ac(x))
    return(tsf)
  } 
  j <- ts(j, start=c(1999, 1), end=c(2018, 12), frequency=12)
      ts.metrics <- ts.features(j)  
    names(ts.metrics)[34] <- "facf"  
  return(ts.metrics)
} 

################################################
# Read lai monthlies (1999-2017) and covariates
# (1999-2010) pre = length(1999:2010)*12 = 1:144 
# (2011-2017) post = post = length(2011:2017)*12 = 145:228
lai <- stack(file.path(ddir,"lai_monthlies.tif"))
  d <- date.seq("1999/01/01", "2017/12/31", step="month")
  d.names <- stringr::str_replace_all(paste0("d", d), "-", ".")
  
settlements <- raster(file.path(getwd(), "data", "dsettlements.tif"))
rds <- raster(file.path(getwd(), "data", "droads.tif"))

################################################
# Create empty model outputs
# Model effect sizes and error rasters
es.raster <- lai[[1]]
  es.raster[] <- rep(NA, ncell(es.raster))  
    writeRaster(es.raster, file.path(getwd(), "acf_results", "effect_sizes.tif"), 
                overwrite=TRUE, options="COMPRESS=LZW")

error.raster <- lai[[1]]
  error.raster[] <- rep(NA, ncell(error.raster))  
    writeRaster(error.raster, file.path(getwd(), "acf_results", "error.tif"), 
                overwrite=TRUE, options="COMPRESS=LZW")

write(paste0("intervention  ", "  lift  ", "  median error  ", "  CATE  ", 
             "  CATT  ", "  CATC"  ), 
	  file.path(getwd(), "acf_results", "cf_results.txt"), 
	  append = FALSE)

################################################
# Create control/experimental units
# Control-case raster classes
#  [1] Social Forestry Areas                [2] KPHP Berau Barat                    
#  [3] Wehea-Kelay                          [4] Padang Luwai Nature Reserve         
#  [5] Muara Kaman Reserve (Nature Reserve) [6] bay reserve                         
#  [7] gulf nature reserve                  [8] Kutai National Park                 
#  [9] Bukit Soeharto Forest Park           [10] Lati Pentangis Forest Park          
# [11] plantation      

# Read experimental units and remove a few problem children
experimental.units <- c("Protected Areas", "IND002.3", "IND002.4", "IND002.5")

eu <- as(st_read(file.path(ddir, "boundaries", "experimental_units.shp")), "Spatial")
  eu <- explode(eu, sp=TRUE)
  # eu <- eu[eu$PROJECT_ID == experimental.units[1] | eu$PROJECT_ID == experimental.units[2],]
    eu$ID <- row.names(eu)

# read plantations
plantations <- as(st_read(file.path(ddir, "boundaries", "concessions.shp")), "Spatial")
  plantations <- explode(plantations, sp=TRUE)
    proj4string(plantations) <- proj4string(eu)

dat <- eu@data[1,]
  dat[1,] <- c("plantation","plantation","2020","999")
    plantations@data <- do.call("rbind", replicate(nrow(plantations), 
	                            dat, simplify = FALSE))

control.units <- rbind(eu, plantations)
  control.units$pid <- as.numeric(factor(control.units$NAME)) 

#### create rasterize control units and create observations
cur <- mask(rasterize(control.units, lai[[1]], field="pid", 
            background=0), lai[[1]]) 
    writeRaster(cur, file.path(getwd(), "acf_results", "eu_raster.tif"), 
                overwrite=TRUE, options="COMPRESS=LZW")

# includes plantations
controls <- as(lai, "SpatialPointsDataFrame")
    names(controls@data) <- d.names
  controls@data <- data.frame(EUID=extract(cur, controls), controls@data)
    controls <- controls[which(controls$EUID == 0 | controls$EUID == 9),] 
      controls@data <- data.frame(EUID=controls@data[,1], treatment=0, 
                                  extract(stack(settlements, rds), controls, cellnumbers=TRUE), 
                                  coordinates(controls), controls@data[,2:ncol(controls)])

na.idx <- unique(which(is.na(controls@data[,4:7]), arr.ind = TRUE)[,1])  
  if(length(na.idx) > 0) {
    controls <- controls[-na.idx,]
  }
 
# plot(cur)
#   plot(controls, pch=20, cex=0.25, add=TRUE) 
 
save.image(file.path(getwd(), "acf_results", "ProcessedData.RData"))
  
#**************************************
#**************************************
# Start of for loop for processing 
# each intervention type
model.data <- list() 
  for(j in unique(eu$PROJECT_ID)) {
    cat("**** Running model for", j, "**** ", "\n")
#**************************************
#**************************************

  # define specific intervention
  # j=unique(eu$PROJECT_ID)[2] 

################################################
# Create interventions and model design matrix

cat("**** Extracting time-series LAI for intervention", j, "\n")
# Subset defined intervention 
  intervention <- exact_extract(lai, as(eu[eu$PROJECT_ID == j,], "sf"), 
                                include_cell = TRUE, include_xy = TRUE)
    intervention <- do.call("rbind", intervention)
      coordinates(intervention) <- ~x+y	  
        intervention@data <- data.frame(EUID=j, treatment=1, cells=intervention@data[,229],
		                                raster::extract(stack(settlements, rds), intervention),
		                                coordinates(intervention),  
										intervention@data[,1:228]) 
          names(intervention@data)[8:ncol(intervention)] <- d.names
		  
na.idx <- unique(which(is.na(intervention@data[,4:7]), arr.ind = TRUE)[,1]) 
  if(length(na.idx) > 0) {
    intervention <- intervention[-na.idx,]
  }

################################################
# kNN and pre-intervention DTW matching
# pre (1999-2010), post (2011-2017) 
intervention.date = "2010-12-31" 
coffset = 7 # column offset of non-temporal variables

pre.idx <- which( d <= intervention.date) + coffset
post.idx <- which( d > intervention.date) + coffset

p = 100     # number of candidate controls, selected using kNN 
nsamp = 3   # number of counterfactuals, selected using DTW  

cat("**** Finding control matches for", j, "\n")

nidx <- RANN::nn2(as.matrix(controls@data[,4:7]),
                  as.matrix(intervention@data[,4:7]), 
			      searchtype="radius", radius = 100000, 
				  k=p)$nn.idx
  nidx <- nidx[-1,]

# Select counterfactual controls using DTW 
control.matches <- data.frame(IDS=0, cells=0, d=0)			
  for(i in 1:nrow(nidx)) {
    cat(i, "in", nrow(nidx), "\n")
    sample.controls <- controls@data[nidx[i,],][,pre.idx]
	ids <- data.frame(IDS=1:nrow(sample.controls), cells=controls@data[nidx[i,],]$cells) 
    reference <- as.numeric(intervention@data[i,][,pre.idx])
    ref.dist <- data.frame(ids, d=unlist(lapply(1:nrow(sample.controls), 
                  FUN=function(j) { dtw(as.numeric(sample.controls[j,]), reference)$distance } ))) 		 
      zero.idx <- which(ref.dist$d <= 0)
  	    if(length(zero.idx) > 0) {
  	      ref.dist <- ref.dist[-zero.idx,]
  	    }
	  control.matches <- rbind(control.matches, ref.dist[which(ref.dist$d %in% 
	                           nth.values(ref.dist$d, N=nsamp, smallest = TRUE)),] )  
  }          
control.matches <- control.matches[-1,]
  control.matches <- control.matches[-which(abs(outliers(control.matches$d)) > 9.9),]
    control.matches <- controls[which(controls$cells %in% unique(control.matches$cells)),]

################################################
# Time-series Feature Engineering 
# Use feature engineering to produce covariates
#   representing characteristics of the timeseries
# https://cran.r-project.org/web/packages/tsfeatures/vignettes/tsfeatures.html
# pre lai value (1999) 1:12 (+7)
# post lai value (2017) 217:228 (+7)

cat("**** Time-series feature Engineering", j, "\n")

coffset = 7  # column offset of non-temporal variables
pre.idx <- which( year(d) == "1999") + coffset     # first year
post.idx <- which( year(d) == "2017") + coffset    # last year

mdata <- rbind(intervention@data, control.matches@data)
  pre.lai = apply(mdata[,pre.idx], MARGIN=1, FUN=max)
    post.lai = apply(mdata[,post.idx], MARGIN=1, FUN=max)
    diff.lai <- pre.lai - post.lai 
  max.lai <- apply(mdata[,8:ncol(mdata)], MARGIN=1, FUN=max) 
  mdata <- data.frame(mdata[,1:coffset], max.lai=max.lai, 
                      pre.lai=pre.lai, post.lai=post.lai, 
					  diff.lai=diff.lai, mdata[,8:ncol(mdata)]) 

mdata <- data.frame(mdata[,1:11], as.data.frame(t(apply(mdata[,12:ncol(mdata)], 
                    MARGIN=1, FUN=process.ts))))
  mdata <- na.omit(mdata)
    mdata$EUID = j
 
################################################
################################################
# Adversarial Causal Forests
# Estimate heterogeneous treatment effects as the
# difference between the expected outcome in the 
# treatment minus the expected outcome in the control
# SMD = ( estimated - control ) / pooled(SD)

# With continuous data the lift is the mean in the treatment minus 
# the mean in the control so, if the treatment yielded a 5.0 
# average LAI and the control yielded 3.5, then the lift would 
# be [5.0 - 3.50 = 1.5]

# For dichotomous (e.g., treated, non-treated), the lift is the 
# probability of the desired outcome in the treatment minus the 
# probability of the desired outcome in the control. 
# For example, if 55% of the treatment showed increased LAI, 
# while 50% in the control showed increased LAI, 
# our lift would be [0.55 - 0.50 = 0.05] 
  
################################################
################################################
# Causal Inference Ensemble model 
# Heterogeneous treatment effect estimation
# Adversarial method using Fast Gradient Sign Method (FGSM)

cat("**** Adversarial Causal Ensemble Model for", j, "\n")
# LAI effect size regression 
( cf <- causal_forest(
  X = model.matrix(~ ., data = mdata[,c(9:ncol(mdata))]),
  Y = mdata$max.lai,
  W = as.numeric(mdata$treatment),
  num.trees = 5001,
  seed = 1839
) )

cat("**** Writing Model Results for", j, "\n")

# generate predictions
preds <- predict( object = cf, estimate.variance = TRUE)
  preds$lower <- preds$predictions - 1.96 * sqrt(preds$variance.estimates)
  preds$upper <- preds$predictions + 1.96 * sqrt(preds$variance.estimates)
  preds$error <- preds$debiased.error + preds$excess.error

# Write model results:
# lift - using upper 50% estimate percentile
# median.error - using sum of: 
#   debiased.error - estimates of the R-loss (Nie and Wager 2017). 
#   excess.error - jackknife estimates of the Monte-carlo error (Wager et al., 2014) 
#                  representing how unstable estimates are.
# CATE - conditional average treatment effect on the full sample
# CATT - conditional average treatment effect on the treated sample
# CATC - conditional average treatment effect on the control sample
lift <- lm(max.lai ~ treatment, mdata[preds$predictions > 
        median(preds$predictions),])
write(paste(j, coef(lift)[2], median(preds$error), 
      average_treatment_effect(cf, target.sample = "all")[1],
	  average_treatment_effect(cf, target.sample = "treated")[1],
      average_treatment_effect(cf, target.sample = "control")[1],
	  sep="  "), 
	  file.path(getwd(), "acf_results", "cf_results.txt"), 
	  append = TRUE)

mdata$predictions  <- preds$predictions
  write.csv(mdata, row.names = FALSE, file.path(getwd(),  
            "acf_results", paste0("model_data_", j, ".csv"))) 

model.data[[j]] <- mdata

# Variable importance
VI <- cf %>% 
  variable_importance() %>% 
    as.data.frame() %>% 
      mutate(variable = colnames(cf$X.orig)) %>% 
        arrange(desc(V1))
VI <- VI[-nrow(VI),]
  VI[,1] <- VI[,1] / max(VI[,1])

################################################
# Assign to prediction raster
pred.idx <- which(mdata$treatment == 1)

es.raster <- raster(file.path(getwd(), "acf_results", "effect_sizes.tif"))
  es.raster[mdata[pred.idx,]$cells] <- preds[pred.idx,]$predictions
    writeRaster(es.raster, file.path(getwd(), "acf_results", "effect_sizes.tif"), 
                overwrite=TRUE, options="COMPRESS=LZW")

error.raster <- raster(file.path(getwd(), "acf_results", "error.tif"))
  error.raster[mdata[pred.idx,]$cells] <- preds[pred.idx,]$error
    writeRaster(error.raster, file.path(getwd(), "acf_results", "error.tif"), 
                overwrite=TRUE, options="COMPRESS=LZW")

#plot(eu[eu$PROJECT_ID == j,])
#  plot(es.raster, add=TRUE)

################################################
# Plot effect sizes with confidence intervals  
y.lim <- c(min(preds$predictions) - min(preds$variance.estimates),
           max(preds$predictions) + max(preds$variance.estimates))
p <- ggplot(data = preds, aes(x = rank(predictions), y = predictions)) + 
            geom_point() + geom_line() +
            labs(x = "Average sample rank", y = "Estimated Treatment Effect")			
  p <- p + geom_ribbon(data = preds, aes(ymin = lower, ymax = upper), 
                       linetype=2, alpha=0.4) + ylim(y.lim)

# Plot Predicted Treatment Effects of 6 top variables
vip <- ggplot(VI, aes(reorder(variable, V1), V1)) +
         geom_bar(stat = "identity") +
           coord_flip()
p1 <- ggplot(mdata, aes_string(x = as.character(VI[1,][2]), y = "predictions")) +
  geom_point() +
    geom_smooth(method = "loess", span = 1) +
      theme_light()
p2 <- ggplot(mdata, aes_string(x = as.character(VI[2,][2]), y = "predictions")) +
  geom_point() +
    geom_smooth(method = "loess", span = 1) +
      theme_light()
p3 <- ggplot(mdata, aes_string(x = as.character(VI[3,][2]), y = "predictions")) +
  geom_point() +
    geom_smooth(method = "loess", span = 1) +
      theme_light()
p4 <- ggplot(mdata, aes_string(x = as.character(VI[4,][2]), y = "predictions")) +
  geom_point() +
    geom_smooth(method = "loess", span = 1) +
      theme_light()
p5 <- ggplot(mdata, aes_string(x = as.character(VI[5,][2]), y = "predictions")) +
  geom_point() +
    geom_smooth(method = "loess", span = 1) +
      theme_light()
p6 <- ggplot(mdata, aes_string(x = as.character(VI[6,][2]), y = "predictions")) +
  geom_point() +
    geom_smooth(method = "loess", span = 1) +
      theme_light()
	  
pdf(file.path(getwd(), "acf_results", paste0(j, ".pdf")), 
    height=8.5, width=11) 
  print(p)
    plot(density(preds$error), main="Error distribution")
      plot(density(preds$excess.error), main="Excess Error distribution")
	    print(vip)
          print(p1)  
            print(p2)
               print(p3)
             print(p4)
           print(p5)
         print(p6)
dev.off()

save.image(file.path(getwd(), "acf_results", paste0("Model_", j, ".RData")))
  remove(cf, preds, mdata, control.matches, intervention, nidx)
    gc()
#**************************************
#**************************************
# End of model(s) for loop
}
#**************************************
#**************************************

model.data <- do.call("rbind", model.data)
  write.csv(model.data, row.names = FALSE, file.path(getwd(),  
            "acf_results", paste0("model_data", ".csv")))

coordinates(model.data) <- ~x+y
  proj4string(model.data) <- proj4string(lai)
    st_write(as(model.data, "sf"), file.path(getwd(),  
            "acf_results", paste0("model_data", ".shp")))

################################################
################################################
# All model results summary and figures
################################################
################################################

es.summary <- function(x, p=0.0001) {
  neutral <- length(x[x > -p & x < p])
  s <- summary(x)
  ne.idx <- which(x %in% x[x > -p & x < p])
    if(length(ne.idx) > 0) x <- x[-ne.idx]
  return(c( negative = length(x[x < 0]), neutral = neutral,
            positive=length(x[x > 0]), s) )
}

################################################
# Summarize global raster results
			
# Extract effect sizes and summarize			
  es <- vector(mode = "list", length = length(unique(eu$PROJECT_ID)))
    names(es) <- unique(eu$PROJECT_ID)  
  for(j in unique(eu$PROJECT_ID)) {
     es[[j]] <- na.omit(do.call("rbind", exact_extract(es.raster, 
	              as(eu[eu$PROJECT_ID == j,], "sf"))))[,1]
  }
  es <- do.call("rbind", lapply(es, es.summary))

# IND002.3 = Social Forestry Areas   
# IND002.4 = Berau Barat  
# IND002.5 = Wehea-Kelay  
# Protected Areas"
es.bar <- data.frame(interventions = rownames(es),
             t(rbind(es[,1] / (es[,1] + es[,2] + es[,3]) * 100,
             es[,2] / (es[,1] + es[,2] + es[,3]) * 100,
             es[,3] / (es[,1] + es[,2] + es[,3]) * 100)))
    names(es.bar)[2:4] <- c("negative","neutral","positive") 
es.bar <- es.bar %>% 
  tidyr::pivot_longer(!interventions, names_to = "effect", values_to = "percent") %>%
    as.data.frame()

es.bar[es.bar == "IND002.3" ,]$interventions <- "Social Forestry Areas (REDD+)" 
  es.bar[es.bar == "IND002.4" ,]$interventions <- "Community Initiative Berau Barat" 
    es.bar[es.bar == "IND002.5" ,]$interventions <- "Wehea-Kelay Ecosystem Essential Area" 

ggplot(es.bar, aes(fill=effect, y=percent, x=interventions)) + 
         geom_bar(position="stack", stat="identity") +
            scale_fill_manual(values = c("coral4", "cornsilk3", "cyan4") ) +
			   coord_flip() +
			     theme(legend.position="none")
			    

################################################
# Summarize local raster results

effect <- function(x) {		
  c( positive = x[3] / (x[1] + x[2] + x[3]) * 100,
     neutral = x[2] / (x[1] + x[2] + x[3]) * 100,
     negative = x[1] / (x[1] + x[2] + x[3]) * 100) 
}
		
# Extract polygon level effect sizes and summarize			
eues <- exact_extract(es.raster, as(eu, "sf"))
  eues <- lapply(eues, function(x) na.omit(as.numeric(x[,1])))
    eues <- lapply(eues, es.summary)
eues.outcome <- do.call("rbind", lapply(eues, effect))
  eu@data <- data.frame(eu@data, eues.outcome)
    eu <- sp.na.omit(eu)  
st_write(as(eu, "sf"), file.path(getwd(), "acf_results", 
         paste0("ExperimentialUnitResults", ".shp")))  



#plot(rank(preds$predictions), preds$predictions, type="n", pch=20,  
#     cex=0.5, ylim = c(-0.02,0.035), xlab = "sample ranks", ylab = "tau")
#  points(rank(preds$predictions), preds$predictions + 1.96 * sigma.hat, pch=20, cex=0.5, col="red")
#  points(rank(preds$predictions), preds$predictions - 1.96 * sigma.hat, pch=20, cex=0.5, col="red")
#  points(rank(preds$predictions), preds$predictions, pch=20, cex=0.5)
