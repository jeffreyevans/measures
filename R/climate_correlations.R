# Bioclim variable descriptions
#   Temp              (BIO1)  Annual Mean Temperature
#   DiurnalRange      (BIO2)  Mean Diurnal Range (Mean of monthly (max temp - min temp))
#   Isothermality     (BIO3)  Isothermality (BIO2/BIO7) (×100) 
#   TempSeasonality   (BIO4)  Temperature Seasonality (standard deviation ×100) 
#   MaxTempWarmMt     (BIO5)  Max Temperature of Warmest Month
#   MinTempColdMt     (BIO6)  Min Temperature of Coldest Month 
#   TempRange         (BIO7)  Temperature Annual Range (BIO5-BIO6)
#   TempWetQt         (BIO8)  Mean Temperature of Wettest Quarter
#   TempDryQt         (BIO9)  Mean Temperature of Driest Quarter
#   TempWarmQt        (BIO10) Mean Temperature of Warmest Quarter
#   TempColdQt        (BIO11) Mean Temperature of Coldest Quarter
#   Precp             (BIO12) Annual Precipitation
#   PrecpWetMt        (BIO13) Precipitation of Wettest Month
#   PrecpDruMt        (BIO14) Precipitation of Driest Month
#   PrecpSeasonality  (BIO15) Precipitation Seasonality (Coefficient of Variation)
#   PrecpWetQt        (BIO16) Precipitation of Wettest Quarter
#   PrecpDryQt        (BIO17) Precipitation of Driest Quarter
#   PrecpWarmQt       (BIO18) Precipitation of Warmest Quarter
#   PrecpColdQt       (BIO19) Precipitation of Coldest Quarter
#   PrecpRange        (NA)    Precipitation Annual Range (BIO16 - BIO17)
#
# Climate change Variables
#   localExtreme 
#   anomalies    
#   velocity     
#   dVelocity    
#   tmax_rate   
#   precip_rate  
#   precip_trend 
#   tmin_trend   
#   tmean_trend  
#   tmax_trend  
#   cmi_trend   
#
suppressMessages(
  lapply(c("sf", "spatialEco", "terra", "duckdb", "geodata", 
           "classInt", "RColorBrewer", "hexbin", "ks"), require, 
		   character.only = TRUE))
sf_use_s2(FALSE)
options(warn=-1)

( cty <- c("bhutan", "costa_rica", "canada", "colombia", "peru", "brazil")[6] )

metric = c("LAI", "FCOV")[1]
clim.type = c("change", "current")[1]

root = file.path("C:/evans/PFP/results", cty)
setwd(root)
  out.dir = file.path(root, "climate")
  dat.dir = file.path("C:/evans/PFP", cty, "data")
  mdl.dir = file.path("C:/evans/PFP", cty, "model", metric)
  clim.dir = file.path("C:/evans/PFP", cty, "climate")
dir.create(out.dir, showWarnings = FALSE)

#********************************************************
# set breaks and colors
clrs <- c("darkred","darkorange1","azure3","aquamarine3","darkgreen")
bks <- c( seq(-1,-0.005, 0.10), c(-0.01, 0.01), seq(0.10, 1, 0.1)) 
  bks.cat <- cut(bks, bks)[-1]	 
    mcls <- colorRampPalette(clrs)(nlevels(bks.cat))
	hex_clrs = colorRampPalette(rev(brewer.pal(11,'Spectral')))
  	  labs <- levels(bks.cat)
        labs <- gsub("[(]", "", labs)
        labs <- gsub("]", "", labs)
		labs <- gsub(",", " to ", labs)

# Bin size for hexbin plots
switch(cty, "brazil" = { bins = 150 }, "bhutan" = { bins = 50 },  
       "canada" = { bins = 50 }, "colombia" = { bins = 100 },  
       "costa_rica" = { bins = 50 }, "peru" = { bins = 100 }) 	

#********************************************************
# Functions
nlcor <- function(x, y, pv = 0.05) {
  g <- mgcv::gam(y ~ s(x))
    g.summ <- summary(g)
  if (g.summ$p.table[4] > pv) {
    nl.corr = 0
  } else {
    nl.corr <- (g$null.deviance - g$deviance)/g$null.deviance
  }
  return(nl.corr)
}

flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  = (cormat)[ut],
    p = pmat[ut]
    )
}

# Xi correlation function 
xi.cor <- function(X, Y, ties = FALSE){
  n <- length(X)
  r <- rank(Y[order(X)], ties.method = "random")
 set.seed(42)
  if(ties){
    l <- rank(Y[order(X)], ties.method = "max")
    return( 1 - n * sum( abs(r[-1] - r[-n]) ) / (2*sum(l*(n - l))) )
  } else {
    return( 1 - 3 * sum( abs(r[-1] - r[-n]) ) / (n^2 - 1) )
  }
}

#********************************************************
# read model results raster
if(metric == "LAI") {
  effect_size <- rast(file.path(root, "lai_effect_sizes.tif"))
  trend <- rast(file.path(root, "lai_tau.tif"))
} else if(metric == "FCOV") {
  effect_size <- rast(file.path(root, "fcov_effect_sizes.tif"))
  trend <- rast(file.path(root, "fcov_tau.tif"))
}

#********************************************************
# read vector data and model results points (cell centers)
rm.vars <- c("DistCrop", "DistDev", "DistWater", "DistStrm", "DistRds", "DistTown", "entropy", "stability",           
             "lumpiness", "hurst", "max_level_shift", "time_level_shift", "max_var_shift", "time_var_shift",
             "max_kl_shift", "time_kl_shift", "crossing_points", "flat_spots", "trend", "spike", "curvature",
             "e_acf1", "e_acf10", "x_acf1", "x_acf10", "diff1_acf1", "diff1_acf10", "diff2_acf1", "diff2_acf10",         
             "seas_acf1", "x_pacf5", "diff1x_pacf5", "diff2x_pacf5", "seas_pacf", "ARCH.LM", "arch_acf", "garch_acf",
             "arch_r2", "garch_r2",  "facf") 

bdy <- st_read(file.path(dat.dir, paste0(cty, ".gpkg")), "boundary") 

db.name <- list.files(mdl.dir, "duckdb$", full.names=TRUE)
  if(length(db.name) > 1) db.name <- grep("results", db.name, value = TRUE)
con <- dbConnect(duckdb::duckdb(), db.name, read_only=TRUE)
  db.table <- dbListTables(con)[1]
  dat <- dbGetQuery(con, paste0("SELECT * FROM ", db.table)) 
dbDisconnect(con, shutdown = TRUE)

dat <- st_as_sf(dat, coords = c("X", "Y"), crs = st_crs(bdy), agr = "constant")
  dat <- dat[,-which(names(dat) %in% rm.vars)]
    e <- extract(c(effect_size, trend), dat)[,-1]
      names(e)[1:3] <- paste0("es_", names(e)[1:3])
        dat <- cbind(dat, e)
          treatment <- dat[dat$NAME != "control",] 
remove(dat)

#********************************************************
# read climate data and assign to effect size points
f <- list.files(clim.dir, full.names = TRUE) 
  switch(clim.type, 
    "change" = { clim <- rast(grep(paste(c("metrics", "trends", "timing"), collapse="|"), f, value = TRUE)) }, 
    "current" = { clim <- rast(grep("bioclime", f, value = TRUE)) }) 	

e <- extract(clim, treatment)[,-1]
  treatment <- cbind(treatment, e)
    treatment <- st_drop_geometry(treatment)

#********************************************************
# Scatter plots
xvars <- names(e)
  rm.idx <- grep("shift", xvars)
    if(length(rm.idx) > 0) xvars <- xvars[-rm.idx]
yvars <- c("es_cpct", "es_current", "es_tsa", "tau_pre", "tau_post")
n = length(xvars) * length(yvars)

ect = 0
for(elu in unique(treatment$elu)) {
#for(elu in unique(treatment$elu)[-c(1:11)]) {
  ect = ect + 1
  cors <- data.frame(matrix(vector(), n, 7)) 
    names(cors) <- c("var1","var2","cor","p","LCI","UCI","xi")
  dir.create(file.path(out.dir, elu), showWarnings = FALSE)
    tdat <- treatment[treatment$elu == elu,]
      if(nrow(tdat) < 1) next 
  pdf(file.path(file.path(out.dir, elu), paste0(metric, "_climate_plots.pdf")), height=10, width=10)
    ct = 0
    for(j in yvars) {
      for(i in xvars) {
  	  ct=ct+1
  	  cat("\n")
	  cat("****  Deriving climate metrics for", elu, "-", ect, "of", length(unique(unique(treatment$elu))),"\n")
  	    flush.console(); Sys.sleep(0.01)
          xy <- cbind(tdat[,j], tdat[,i])
  		  colnames(xy) <- c(j,i)
            na.idx <- as.data.frame(which(is.na(xy), arr.ind = TRUE))[,1]
		      if(length(na.idx) > 0) xy <- xy[-na.idx,] 
        if(nrow(xy) < 1) next
        if(length(unique(xy[,1])) / nrow(xy) <= 0.10 | length(unique(xy[,2])) / nrow(xy) <= 0.10) next
  	  cat("****  Creating HexBin plots for", j, "and", i, "across all effect sizes", "\n")
	  cat("****    ", ct, "of", n, "\n")
  	    flush.console(); Sys.sleep(0.01)
  	  bin <- hsmooth(hexbin(xy[,j], xy[,i], xbins = bins, ID=TRUE), wts = c(24,12,6))
          idx <- which(bin@count < quantile(bin@count,p=0.50))
          bin@xcm <- bin@xcm[-idx]
          bin@ycm <- bin@ycm[-idx]
          bin@count <- bin@count[-idx]
          bin@cell <- bin@cell[-idx]
  	  cat("****  Deriving KDE and dropping low-density regions", "\n")                
        fhat <- suppressWarnings(ks::kde(xy, h=Hscv(na.omit(xy))))
          p = contourLevels(fhat, prob=c(0.75))
            est <- suppressWarnings(ks::kde(xy, h=Hscv(na.omit(xy)), eval.points = xy[,1:2])$estimate)
          	  if(length(na.idx) > 0) est <- insert.values(est, NA, na.idx)
  			    if(round(p,5) > 0) { 
  			      xy <- as.data.frame(xy)[which(est > p),]
  			      if(nrow(xy) == 0) xy <- cbind(tdat[,j], tdat[,i]) 
  			  }
  	    cat("****  Fitting GLM for", j, "and", i, "\n")
  		cat("\n")
          flush.console(); Sys.sleep(0.01)
        if(nrow(xy) < 10) next				  
        if(length(unique(xy[,1])) / nrow(xy) <= 0.10 | length(unique(xy[,2])) / nrow(xy) <= 0.10) next		  
		g <- tryCatch( expr = { mgcv::gam(xy[,j] ~ s(xy[,i])) }, error = { next })
		  if(round(summary(g)$dev.expl,4)*100 < 5) next 		 
        gplot.hexbin(bin, colramp = hex_clrs, legend = FALSE, 
                     main = paste0("Scatter-plot density of ", j, " and ", i, " for all effects"), 
                     xlab = j, ylab = i)		  
  		plot(g, main = paste0("Non-linear fit of ", j, " and ", i, " for all effects"),
  		     xlab = i, ylab = j, sub=paste0("R-square = ", round(summary(g)$r.sq,3),
               " Deviance explained = ", round(summary(g)$dev.expl,4)*100,"%"), rug=FALSE)		   
          # Correlations
          s <- round(unlist(cor.test(xy[,j], xy[,i], type="spearman")[c("estimate","p.value","conf.int")]), 6)
            cors[ct,][1:2] <- c(j,i)
  		      cors[ct,][3:7] <- c(s, xi.cor(xy[,j], xy[,i]))
      }
    }
  dev.off()
  write.csv(cors, file.path(file.path(out.dir, elu),  paste0(metric, "_climate_correlations.csv")))
}
