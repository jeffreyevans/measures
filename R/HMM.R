# Hidden Markov Model for forest disturbance
load("C:/evans/measures/Indonesia/CIA_model_data.RData")

library(depmixS4)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(sp)
library(raster)
library(terra)
library(spatialEco)
library(rgdal)
library(rts)
library(dtw)

setwd("C:/evans/measures/Indonesia")

eu <- readOGR(file.path(getwd(), "boundaries"), "experimental_units")
  #eu.sub <- eu[eu$NAME == "Wehea-Kelay",]
  #eu.sub <- as(extent(12942395, 13007403, 70520.73, 121082.1), "SpatialPolygons") 

pre.intervention <- which(lai_dates > as.Date("1998-12-31") &
                          lai_dates <= as.Date("2007-12-31"))
post.intervention <- which(lai_dates >= as.Date("2009-01-01") &
                           lai_dates < as.Date("2018-01-10"))
sub.dates <- post.intervention

lai.sub <- rast(lai[[sub.dates]])
  lai.sub <- mask(crop(lai.sub, ext(eu.sub)),vect(eu.sub))
  #dat <- as.numeric(terra::extract(lai.sub, cbind(13000829, 110061.2)))[-1]
  #dat <- as.numeric(terra::extract(lai.sub, cbind(12976892, 174961.3)))[-1]
  #laiv <- as(stack(lai.sub), "SpatialPixelsDataFrame")
  #  laiv <- na.omit(laiv@data)

#********************************************** 
# Build expected distribution
# 1. Calculate proportion of temporal observations lai >= 6
# 2. Mask out observations < 0.80
# 3. Calculate min, max and median column-wise distributions
# 4. Calculate trend of median using Seasonal detrending using 
#    Bayesian Estimator of Abrupt change, Seasonal change, 
#    and Trend (BEAST algorithm)
lai.f <- app(lai.sub, function(x, p=6) {
                x.sub <- na.omit(x[x > p])
                return(c( length(x.sub) / length(x) )) 
			  })
	lai.f[lai.f <= 0.8] <- NA		  
lai.d <- as( stack(mask(lai.sub, lai.f)), "SpatialPixelsDataFrame")
  lai.d <- na.omit(lai.d@data)

dmin <- apply(lai.d, 2, min)  
  dmed <- apply(lai.d, 2, median)
    dmax <- apply(lai.d, 2, max)
      expected <- Rbeast::beast(dmed, period=36)$t

# Plot bounds and expected based on trend  
plot(lai_dates, dmed, type="l", 
     ylim=range(c(dmin,dmax)))
  lines(lai_dates, expected, col="red")
    lines(lai_dates, dmin, col="grey" ,lty=3)  
      lines(lai_dates, dmax, col="grey",lty=3)  
   
#**********************************************    
# Draw a random obs to test
dat <- as.numeric(sampleRandom(stack(lai), 1))
  dat <- data.frame(ID=1:length(dat), 
           date=lai_dates, y=dat) 
  x <- dat$y
  
#************************************************* 
# Hidden State Markov Model
# predict the states by estimating the posterior
# x        Numeric time series vector
# s        Number of states
# min.run  Minimum number of sequential observations 
#          in a state run. Is predicated on time-step 
#          of time-series ie., if monthly than min.run=12 
#          would represent a minimum run of 12 months
mmhs <- function(x, s = 2, min.run = 6) {
  fit.mod <- fit(depmixS4::depmix(x ~ 1, nstates = s, 
                 family = gaussian(), ntimes=length(x)))
  e <- depmixS4::posterior(fit.mod)
    e$run <- sequence(rle(as.character(e$state))$lengths)
    state.class <- rep(NA,length(e$state))
      state.class[1] <- 1
      j=1
      for(i in 2:(length(e$state)-1)) {
        if(e$state[i] == e$state[i+1]) {
         state.class[i] <- j
       } else {
         j=j+1
         state.class[i] <- j
       }
      }
	  state.class[length(state.class)] <- state.class[length(state.class)-1]
        period.lengths <- rle(as.character(e$state))$lengths
          idx <- which(period.lengths < min.run)		
  if(length(idx) > 0) {
    period.lengths[idx] <- NA
      state.class[which(state.class %in% idx)] <- NA
	    na.idx <- which(is.na(state.class))
    state.class[na.idx] <- ceiling(approx(x = state.class, xout = na.idx, 
                                   method = "constant")$y)
	  period.lengths <- rle(as.character(state.class))$lengths							   
  }	
  n.periods <- length(unique(period.lengths[!is.na(period.lengths)]))
  return( list(idx = 1:length(x), y = x,
               duration = mean(period.lengths, na.rm=TRUE),
               number.periods = n.periods,
  			   classes = state.class) )
  } 
states <- mmhs(x) 


# Calculate seasonal de-trended 
trend <- Rbeast::beast(x, period=36)$t

# create color vector
nc = length(unique(states$classes))
  classes <- unique(states$classes)
    cls <- topo.colors(10)
       cls <- sample(cls, nc, replace=TRUE)
my.clr <- states$classes
  for(i in 1:nc) { my.clr[my.clr == classes[i]] <- cls[i] }
  par(mfrow=c(2,1))
    plot(lai_dates, x, type="l", xlab="", ylab="raw")
      points(lai_dates, x, col=my.clr, pch=20)	
    plot(lai_dates, expected, type="l", xlab="", ylab="trend") 
	

    plot(lai_dates, expected, type="l", xlab="", ylab="Expected distribution",
	     ylim=range(c(x,expected))) 
      lines(lai_dates, x, lty=3)
	    points(lai_dates, x, col=my.clr, pch=20, cex=0.75)
	
  par(mfrow=c(2,1))
    plot(dat$date, dat$trend, type="l", xlab="", ylab="trend") 
      points(dat$date, dat$trend, col=my.clr, pch=20)
    plot(dat$date, dat$S1,type="l", col="red", lwd=2, ylab="state probs",
         xlab="")
      lines(dat$date, dat$S2, col="blue", lwd=2)  
  par(mfrow=c(1,1))   
    plot(d,type="two",off=1, ylab="DTW", xlab="")	  
dev.off()  


#***********************************************
# Format data  
dat <- data.frame(dat, est.states[,4:ncol(est.states)])
     period.lengths <- rle(as.character(dat$state))$lengths
        n.periods <- length(period.lengths)
         
sdat <- split(dat, factor(est.states$class)) 

# Creates equal length vectors based on adding NA's 
# then imputing the missing values 
equal.length <- function(x, n) {
  impute.loess <- function(y, s = 0.2) {
    x.length = length(y)
    if(length(y[!is.na(y)]) < 6) {
      warning("Fewer than 6 real-value observations, assigning NA")
        y <- rep(NA, x.length)
    } else {
      x <- 1:x.length
        p <- suppressWarnings(stats::loess(y ~ x, span = s, 
                              data.frame(x = x, y = y)))
        na.idx <- which(is.na(y))
          if(length(na.idx) > 1) {
            y[na.idx] <- stats::predict(p, data.frame(x = na.idx))
          }
      }
    return(y)
  }
  if(length(x) == n)
    stop("Vectors is already = ", n, " length")
	  d = n - length(x)
	  x <- spatialEco::insert.values(x, NA, sample(3:n,d))
    x <- impute.loess(x)
  return(x)   
}

n = max(unlist(lapply(sdat, nrow)))
trends <- list()
  for(i in 1:length(sdat)) {
    if(nrow(sdat[[i]]) == n) {
	  trends[[i]] <- as.numeric(sdat[[i]]$trend)
	} else {
	  d = n - nrow(sdat[[i]])
      trends[[i]] <- equal.length(sdat[[i]]$trend[-1], n)
    }
  }
xy <- do.call(cbind, trends)

#****************************************************
# Dynamic Time Warping distances
u <- combn(c(1:ncol(xy)), 2)
dtw.dist <- vector()
  for(i in 1:ncol(u)) {
    dtw.dist[i] <- dtw::dtw(xy[,u[,i][1]], xy[,u[,i][2]], 
                            window.type = "sakoechiba",
    		                window.size = 1, 
							distance.only=TRUE)$distance / 
				            abs(sum(xy[,u[,i][1]]))
  }
  
# dtw example with plot
x <- trends[[1]]
y <- trends[[2]]
  na.idx <- c(which(is.na(x)),which(is.na(x)))
    if(length(na.idx) > 0) {
	  x <- x[-na.idx]
	  y <- y[-na.idx]
	}
d <- dtw::dtw(x[-194], y[-194], window.type = "sakoechiba", 
			  keep=TRUE, step=asymmetric,
              open.end=TRUE,open.begin=TRUE,
		      window.size = 1)
  #plot(d, type="two", off=1)
  
#*********************************************
# Plot results

# create color vector
nc = length(unique(dat$class))
  classes <- unique(dat$class)
    cls <- topo.colors(20)
       cls <- sample(cls,nc)
my.clr <- dat$class 
  for(i in 1:nc) { my.clr[my.clr == classes[i]] <- cls[i] }

  par(mfrow=c(2,1))
    plot(dat$date, dat$y, type="l", xlab="", ylab="raw")
    plot(dat$date, dat$trend, type="l", xlab="", ylab="trend")
  par(mfrow=c(2,1))
    plot(dat$date, dat$trend, type="l", xlab="", ylab="trend") 
      points(dat$date, dat$trend, col=my.clr, pch=20)
    plot(dat$date, dat$S1,type="l", col="red", lwd=2, ylab="state probs",
         xlab="")
      lines(dat$date, dat$S2, col="blue", lwd=2)  
  par(mfrow=c(1,1))   
    plot(d,type="two",off=1, ylab="DTW", xlab="")	  
dev.off()  
  
  
#***********************************************
# calculate resilience metrics    

tau <- vector() 
slp <- vector() 
  for(i in 1:length(sdat)) {
    tau[i] <- kendall(sdat[[i]]$trend)[2]
	slp[i] <- kendall(sdat[[i]]$trend)[1]
  }

# Pairwise Kullback-Leibler divergence
u <- combn(c(1:ncol(xy)), 2)
kl <- c(NA)
  for(i in 1:ncol(u)) {
    kl[i] <- kl.divergence(xy[,c(u[,i][1],u[,i][2])])[3:3]
  }  
#plot(sort(kl), type="b")

# Pairwise Hausdorff Distance
hasd <- vector()
  for(i in 1:ncol(u)) {
    hasd[i] <- pracma::hausdorff_dist(xy[,u[,i][1]],xy[,u[,i][2]])
  }
#plot(sort(hasd), type="b")

# cp.idx <- c(changepoint.np::cpt.np(dat$trend, minseglen=10)@cpts)
# par(mfrow=c(2,1))
#   plot(dat$y, type="l")
#     points(cp.idx, dat$y[cp.idx], pch=19)
#   plot(dat$trend, type="l")
#     points(cp.idx, dat$trend[cp.idx], pch=19)  

#****************************************************
# Dynamic Time Warping distances
dtw.dist <- vector()
  for(i in 1:ncol(u)) {
    dtw.dist[i] <- dtw::dtw(xy[,u[,i][1]], xy[,u[,i][2]], 
                            window.type = "sakoechiba",
    		                window.size = 1, 
							distance.only=TRUE)$distance / 
				            abs(sum(xy[,u[,i][1]]))
  }

#*********************************************
# Plot results ggplot
ns.lab <- as.character(1:length(grep("S", names(est.states))))
mycols <- c("darkmagenta", "blue")
g0 <- (ggplot(dat, aes(x = date, y = y)) + geom_line() +
       theme(axis.ticks = element_blank(), 
	   axis.title.y = element_blank()) +
	   labs(y = "Raw time series")) %>% ggplotGrob	
g1 <- (ggplot(dat, aes(x = date, y = trend)) + geom_line() +
       theme(axis.ticks = element_blank(), 
	   axis.title.y = element_blank()) +
	   labs(y = "De-trended time series")) %>% ggplotGrob	
g2 <- (ggplot(est.states, aes(x = ID, y = state),
       fill = ns.lab, col = ns.lab) +
       geom_bar(stat = "identity") +
       theme(axis.ticks = element_blank(), axis.text.y = element_blank()) +
       labs(y = "Estimated State")) %>% ggplotGrob	   
g3 <- (ggplot(est.long, aes(x = ID, y = value, col = variable)) + geom_line() +
       theme(axis.ticks = element_blank(), axis.text.y = element_blank()) + 
       labs(y = "Posterior Prob.")) %>% ggplotGrob()	      	   
	g1$widths <- g2$widths
grid.arrange(g0, g1, g2, g3, widths = 1, nrow = 4)  
  
  