##################################################
##################################################
#             Causal Impact Model 
#     Copernicus 1000m resolution LAI data 
##################################################
##################################################
# [1] "Protected Areas",
# [2] "IND002.3_Berau_Community_Initative_Social_Forestry_Areas",
# [3] "IND002.4_Berau_Community_Initative_Barat_KPH_FMU_Model_Area",
# [4] "IND002.5_Berau_Community_Initative_Wehea_Kelay_Ecosystem_Essential_Area"
	
load("C:/evans/measures/Indonesia/CIA_model_data.RData")
  
#### fetch specific objects in a lazy load
# e = local({load("C:/evans/measures/Indonesia/CIA_model_data.RData"); environment()})
#   tools:::makeLazyLoadDB(e, "CIA_model_data")
#     lai.controls 
#     lai_dates
#     lai.all  
  
library(sp)
library(raster)
library(spatialEco)
library(rgdal)
library(prophet)
library(CausalImpact)
library(dtw)
#library(interventions)
library(MarketMatching)
library(dplyr)
library(lubridate)


setwd("C:/evans/measures/Indonesia")

#### Read experimental units and remove a few problem children
eu <- eu[eu$PROJECT_ID == experimental.units[1] | eu$PROJECT_ID == experimental.units[2],]
  eu$ID <- row.names(eu)

#### Remove controls occurring in oil-palm concessions
palm <- readOGR(file.path(getwd(), "boundaries"), "concessions")
  proj4string(palm) <- proj4string(eu)

o = over(lai.controls, palm)
lai.controls <- lai.controls[-which(!is.na(o$OBJECTID)),]

#  plot(palm,
#       xlim=c(12677311, 12819332),
#       ylim=c(42294.43, 162959.1))
#  plot(lai.controls,pch=20,cex=0.2,add=TRUE)

###############################################
#### Smooth time-series using loess smoothing 
impute.loess <- function(y, x.length = NULL, s = 0.05, sdata = TRUE) {
  if(is.null(x.length)) { x.length = length(y) }
    x <- 1:x.length
	  p <- stats::loess(y ~ x, span = s, data.frame(x=x, y=y))
	if(sdata == TRUE) {
	  y <- stats::predict(p, x)
	} else {
    na.idx <- which( is.na(y) )
      if( length(na.idx) > 1 ) {
        y[na.idx] <- stats::predict(p, data.frame(x=na.idx))
	    }
  }  
  return(y)
}
   
# lai.controls@data <- as.data.frame(t(apply(lai.controls@data, MARGIN=1, FUN=impute.loess)))
# lai.all@data <- as.data.frame(t(apply(lai.all@data, MARGIN=1, FUN=impute.loess)))
     
#### Truncate time-period
idx <- which(format(as.Date(lai_dates, format="%d/%m/%Y"),"%Y") < "2002" |
             format(as.Date(lai_dates, format="%d/%m/%Y"),"%Y") > "2017")
  lai_dates <- lai_dates[-idx]
    lai.controls@data <- lai.controls@data[,-idx]
	  lai.all@data <- lai.all@data[,-idx]
  
###################
# model parameters
###################
n = round(nrow(lai.controls) * 0.010, 0)  # Number of random pixel samples to test with DTW model
idate = "2008-01-01"                      # Intervention date
ncontrol = 100                            # Number of controls

###############################################
#### Extract raster values, calculate summary 
####   statistic for each experimental unit and
####   create a data.frame of each pixel trend
ldf <- point.in.poly(lai.all, eu)
  ldf <- data.frame(ID=ldf$ID, ldf@data[1:length(lai_dates)])
    names(ldf)[2:ncol(ldf)] <- as.character(lai_dates)
	  ldf <- na.omit(ldf)
	  
#### Random subsample
#ldf <- ldf[sample(1:nrow(ldf), nrow(ldf) * 0.10),] 

###############################################
###############################################
#### Start loop to process all pixels in each
####  experimental unit

cia.coefficents <- list()
cia.average <- list()
cia.cummulative <- list()
trend <- list()

#### start model(s) loop
  for(i in 1:nrow(ldf)) {
    lai.dat <- data.frame(ID = "intervention", 
                          date = lai_dates, 
                          lai = as.numeric(ldf[i,][2:ncol(ldf)]))

###############################################
#### Create distance-weighted random sample 
####   for experimental unit associated with
####   given i pixel
euid <- ldf[i,]$ID
d <- rgeos::gDistance( eu[which(euid %in% row.names(eu)),], lai.controls, byid=TRUE)
  p <- ((d - max(d)) * -1) + min(d)
    p <- p / sum(p)
  rs.lai <- data.frame(ID = paste0("control",1:n), 
               lai.controls@data[sample(1:nrow(lai.controls@data),n, prob = p),] )			   
  for(j in 1:n){
    lai.dat <- rbind(lai.dat, data.frame(ID = rep(as.character(rs.lai[j,][1]), 
	                 length(lai_dates)), date = lai_dates,
	                 lai = as.numeric(rs.lai[j,][-1])) )  
  }
    
##################################################
##################################################
# Dynamic Time Warping matching
##################################################
##################################################

cat("Testing", round(nrow(lai.controls) * 0.010, 0), "(1% sample) temporal distributions for",  
    ncontrol, "controls in establishing priors for observation", i, "\n")
mm <- best_matches(data = lai.dat,
                   id_variable = "ID",
                   markets_to_be_matched="intervention",
                   date_variable = "date",
                   matching_variable = "lai",
                   parallel = TRUE,
                   warping_limit = 1, 
                   dtw_emphasis = 1, 
                   matches = ncontrol, 
                   start_match_period = min(lai_dates),
                   end_match_period = idate)
			      
##################################################
##################################################				   
# LAI Causal Impact Bayesian Model using BSTS 
#   (Bayesian structural time-series models) to 
#   account for seasonal periodicity with spike-and-slab 
#   priors 
##################################################
##################################################

trend[[i]] <- data.frame(ID=ldf[i,]$ID, round(kendall(as.numeric(lai.dat[lai.dat$ID == "intervention",]$lai)),5))

#x <- data.frame(lai=as.numeric(lai.dat[lai.dat$ID == "intervention",]$lai),  
#                date=lai_dates) 
#x <- x %>% mutate(month = format(date, "%m"), year = format(date, "%Y")) %>%
#               group_by(month, year) %>% 
#			     summarise(lai = median(lai)) %>% 
#				   as.data.frame()
#x <- as.ts(reshape2::acast(x, year~month, value.var="lai"), 
#        start=c(2002, 1), end=c(2017, 1), frequency=1)

cia <- inference(matched_markets = mm,
                 bsts_modelargs = list(niter = 1000, nseasons = 36), 
                 alpha = 0.05, 
                 prior_level_sd = 0.01,
                 control_matches = 10,
                 analyze_betas = FALSE,
                 test_market = "intervention",
                 end_post_period = lai_dates[length(lai_dates)])
			
cof <- as.data.frame(rbind(as.numeric(cia$Coefficients[,2][1:10])))
  names(cof) <- paste0("coefficient", 1:10)
cia.coefficents[[i]] <- data.frame(ID=ldf[i,]$ID, 
                                   AbsoluteEffect=cia$AbsoluteEffect,
                                   AbsoluteEffectLower=cia$AbsoluteEffectLower,
                                   AbsoluteEffectUpper=cia$AbsoluteEffectUpper,
                                   RelativeEffect=cia$RelativeEffect*100,
					               RelativeEffectLower=cia$RelativeEffectLower*100,
                                   RelativeEffectUpper=cia$RelativeEffectUpper*100,
					               p=cia$TailProb, MAPE=cia$PrePeriodMAPE*100, 
                                   DW=cia$DW, cof)
cia.average[[i]] <- cia$CausalImpactObject$summary[1,] 
cia.cummulative[[i]] <- cia$CausalImpactObject$summary[2,]  
}

cia.coefficents <- do.call(rbind, cia.coefficents)
cia.average <- do.call(rbind, cia.average)
cia.cummulative <- do.call(rbind, cia.cummulative)
trend <- do.call(rbind, trend)


# cia.model <- cia$CausalImpactObject				 			 
#cia$PlotActualVersusExpected
#cia$PlotCumulativeEffect
#cia$PlotPointEffect
#cia$PlotActuals
# summary(cia$CausalImpactObject, "report")  
