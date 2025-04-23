suppressMessages(lapply(c("sf", "spatialEco", "terra", "ggplot2"), 
		         require, character.only = TRUE))

( cty <- c("brazil", "bhutan", "canada", "colombia", "costa_rica", "peru")[1] )

root = file.path("C:/evans/PFP/results", cty)
setwd(root)
  out.dir = file.path("C:/evans/PFP/results", cty)
  dat.dir = file.path("C:/evans/PFP", cty, "data")

#********************************************************
# Lables and colors 
pretty_cut <- function(values,breaks,format=function(d)d,
                       binding="to",spacing=" ",
                       under_text="<",over_text=">",
                       ...){
  labels <- seq(1,length(breaks)-1) %>%
    lapply(function(i){
      a=breaks[i]
      b=breaks[i+1]
      if (is.integer(b)) b=b+1
      if (is.infinite(a) & a<0) {
        text=paste0(under_text,spacing,format(b))
      } else if (is.infinite(b) & b>0) {
        text=paste0(over_text,spacing,format(a))
      } else {
        ta=paste0(format(a),spacing,binding,spacing,format(b))
      }
    }) |>
    unlist()
  cut(values,breaks=breaks,labels=labels,...)
}
clrs <- c("darkred","darkorange1","azure3","aquamarine3","darkgreen")
bks <- c( seq(-1,-0.005, 0.10), c(-0.01, 0.01), seq(0.10, 1, 0.1)) 
  hte.bks <- pretty_cut(bks, bks)[-1]	 
  hte.mcls <- colorRampPalette(clrs)(nlevels(hte.bks))
  hte.labs <- levels(hte.bks)

#********************************************************
# read rasters
lai <- rast("lai_effect_sizes_classified.tif", lyr=c(1,2,5))
  for(i in 1:nlyr(lai)){
    l <- levels(lai[[i]])[[1]]
    levels(lai[[i]]) <- NULL
    v <- unique(lai[[i]])[,1]
    coltab(lai[[i]]) <- data.frame(value = v, col = hte.mcls[v])
    levels(lai[[i]]) <- l
  }  
fcov <- rast("fcov_effect_sizes_classified.tif", lyr=c(1,2,5))
  for(i in 1:nlyr(fcov)){
    l <- levels(fcov[[i]])[[1]]
    levels(fcov[[i]]) <- NULL
    v <- unique(fcov[[i]])[,1]
    coltab(fcov[[i]]) <- data.frame(value = v, col = hte.mcls[v])
    levels(fcov[[i]]) <- l
  }  

lai.trend <- rast("lai_tau_classified.tif", lyr=c(3,4))
  names(lai.trend) <- c("post_tau", "pre_tau")
  for(i in 1:nlyr(lai.trend)){
    l <- levels(lai.trend[[i]])[[1]]
    levels(lai.trend[[i]]) <- NULL
    v <- unique(lai.trend[[i]])[,1]
    coltab(lai.trend[[i]]) <- data.frame(value = v, col = hte.mcls[v])
    levels(lai.trend[[i]]) <- l
  }  

fcov.trend <- rast("fcov_tau_classified.tif", lyr=c(3,4))
  names(fcov.trend) <- c("post_tau", "pre_tau")
  for(i in 1:nlyr(fcov.trend)){
    l <- levels(fcov.trend[[i]])[[1]]
    levels(fcov.trend[[i]]) <- NULL
    v <- unique(fcov.trend[[i]])[,1]
    coltab(fcov.trend[[i]]) <- data.frame(value = v, col = hte.mcls[v])
    levels(fcov.trend[[i]]) <- l
  }  

#********************************************************
# Read boundaries
bdy <- st_read(file.path(dat.dir, paste0(cty, ".gpkg")), "boundary")
if(any(cty %in% c("colombia", "peru", "brazil"))){ 	  
  pa <- st_read(file.path(dat.dir, paste0(cty, ".gpkg")), "protected_areas")
    st_geometry(pa) <- "geometry" 
      if(cty == "brazil") pa <- pa[-which(pa$DESIG %in% "Indigenous Area"),]
  	pa$type <- "protected Area"
  pfp <- st_read(file.path(dat.dir, paste0(cty, ".gpkg")), "PFP")
    st_geometry(pfp) <- "geometry"
  	pfp$IUCN <- "Not Applicable"
	pfp$type <- "PFP"
    interventions <- st_read(file.path(dat.dir, paste0(cty, ".gpkg")), "intervention_boundaries") 
} else {
  pa <- st_read(file.path(dat.dir, paste0(cty, ".gpkg")), "protected_areas")
    st_geometry(pa) <- "geometry"
  interventions <- pa	
}

#********************************************************
# Plot effect size maps

dev.new()

# lai current
plot(lai[["current"]], legend=FALSE, maxcell=50000000, box=NA, axes= FALSE)
  plot(st_geometry(bdy), add=TRUE)
    plot(st_geometry(interventions), add=TRUE, lwd=0.75)

# lai tsa
plot(lai[["tsa"]], legend=FALSE, maxcell=50000000, box=NA, axes= FALSE)
  plot(st_geometry(bdy), add=TRUE)
    plot(st_geometry(interventions), add=TRUE, lwd=0.75)

# lai change
plot(lai[["change"]], legend=FALSE, maxcell=50000000, box=NA, axes= FALSE)
  plot(st_geometry(bdy), add=TRUE)
    plot(st_geometry(interventions), add=TRUE, lwd=0.75)


# fcov current
plot(fcov[["current"]], legend=FALSE, maxcell=50000000, box=NA, axes= FALSE)
  plot(st_geometry(bdy), add=TRUE)
    plot(st_geometry(interventions), add=TRUE, lwd=0.75)

# fcov tsa
plot(fcov[["tsa"]], legend=FALSE, maxcell=50000000, box=NA, axes= FALSE)
  plot(st_geometry(bdy), add=TRUE)
    plot(st_geometry(interventions), add=TRUE, lwd=0.75)

# fcov change
plot(fcov[["change"]], legend=FALSE, maxcell=50000000, box=NA, axes= FALSE)
  plot(st_geometry(bdy), add=TRUE)
    plot(st_geometry(interventions), add=TRUE, lwd=0.75)


#********************************************************
# Effect size bar plots
f <- data.frame(freq(lai[["current"]]))
  f$count <- f$count / sum(f$count)
barplot(f$count, col = hte.mcls, ylim = c(0, (max(f$count)+0.05)), space=0, beside = TRUE)

f <- data.frame(freq(lai[["tsa"]]))
  f$count <- f$count / sum(f$count)
barplot(f$count, col = hte.mcls, ylim = c(0, (max(f$count)+0.05)), space=0, beside = TRUE)

f <- data.frame(freq(lai[["change"]]))
  f$count <- f$count / sum(f$count)
barplot(f$count, col = hte.mcls, ylim = c(0, (max(f$count)+0.05)), space=0, beside = TRUE)


f <- data.frame(freq(fcov[["current"]]))
  f$count <- f$count / sum(f$count)
barplot(f$count, col = hte.mcls, ylim = c(0, (max(f$count)+0.05)), space=0, beside = TRUE)

f <- data.frame(freq(fcov[["tsa"]]))
  f$count <- f$count / sum(f$count)
barplot(f$count, col = hte.mcls, ylim = c(0, (max(f$count)+0.05)), space=0, beside = TRUE)

f <- data.frame(freq(fcov[["change"]]))
  f$count <- f$count / sum(f$count)
barplot(f$count, col = hte.mcls, ylim = c(0, (max(f$count)+0.05)), space=0, beside = TRUE)


#********************************************************
# Plot trend maps

dev.new()
par(mfrow=c(1,2))

# lai current
plot(lai.trend[[2]], legend=FALSE, maxcell=50000000, box=NA, axes= FALSE)
  plot(st_geometry(bdy), add=TRUE)
    plot(st_geometry(interventions), add=TRUE, lwd=0.75)
plot(lai.trend[[1]], legend=FALSE, maxcell=50000000, box=NA, axes= FALSE)
  plot(st_geometry(bdy), add=TRUE)
    plot(st_geometry(interventions), add=TRUE, lwd=0.75)

# fcov current
plot(fcov.trend[[2]], legend=FALSE, maxcell=50000000, box=NA, axes= FALSE)
  plot(st_geometry(bdy), add=TRUE)
    plot(st_geometry(interventions), add=TRUE, lwd=0.75)
# lai tsa
plot(fcov.trend[[1]], legend=FALSE, maxcell=50000000, box=NA, axes= FALSE)
  plot(st_geometry(bdy), add=TRUE)
    plot(st_geometry(interventions), add=TRUE, lwd=0.75)
  

#********************************************************
# Trend bar plots

# LAI
f <- data.frame(freq(lai.trend, usenames=TRUE))
  names(f) <- c("metric", "effect_size", "pct")
    f$metric <- gsub("pa_tau_post", "tau_post", f$metric)
	f$metric <- gsub("pa_tau_pre", "tau_pre", f$metric)
    f$pct[which(f$metric == "tau_post")] <- f[which(f$metric == "tau_post"),]$pct/ sum(f[which(f$metric == "tau_post"),]$pct) 
    f$pct[which(f$metric == "tau_pre")] <- f[which(f$metric == "tau_pre"),]$pct/ sum(f[which(f$metric == "tau_pre"),]$pct) 
f <- reshape(f, idvar = "metric", timevar = "effect_size", direction = "wide")
names(f)[-1] <- gsub("pct.", "es ", names(f)[-1])

barplot(as.numeric(f[1,][-1]), col = "#006400", ylim = c(0, (max(as.numeric(f[1,][-1]))+0.05)), space=0, beside = TRUE)
  barplot(as.numeric(f[2,][-1]),col = adjustcolor("coral4", 0.6), space=0, add = TRUE)
legend("topleft", legend=c("Tau pre-intervention", "Tau post-intervention"),
       fill = c(adjustcolor("coral4", 0.6),"#006400"))

# FCOV
f <- data.frame(freq(fcov.trend, usenames=TRUE))
  names(f) <- c("metric", "effect_size", "pct")
    f$metric <- gsub("pa_tau_post", "tau_post", f$metric)
	f$metric <- gsub("pa_tau_pre", "tau_pre", f$metric)
    f$pct[which(f$metric == "tau_post")] <- f[which(f$metric == "tau_post"),]$pct/ sum(f[which(f$metric == "tau_post"),]$pct) 
    f$pct[which(f$metric == "tau_pre")] <- f[which(f$metric == "tau_pre"),]$pct/ sum(f[which(f$metric == "tau_pre"),]$pct) 
f <- reshape(f, idvar = "metric", timevar = "effect_size", direction = "wide")
names(f)[-1] <- gsub("pct.", "es ", names(f)[-1])

barplot(as.numeric(f[1,][-1]), col = "#006400", ylim = c(0, (max(as.numeric(f[1,][-1]))+0.05)), space=0, beside = TRUE)
  barplot(as.numeric(f[2,][-1]),col = adjustcolor("coral4", 0.6), space=0, add = TRUE)
legend("topleft", legend=c("Tau pre-intervention", "Tau post-intervention"),
       fill = c(adjustcolor("coral4", 0.6),"#006400"))


#********************************************************
# Change in trend table
lai.chg <- rast("lai_tau_change.tif")
  d <- na.omit(data.frame(c(lai, lai.chg)))
    ft <- lapply(c("current", "tsa"), \(i) { 
      f <- data.frame(metric = i, table(d[,"lai_change_type"], d[,i]))
        names(f) <- c("metric",  "trend", "effect_size", "count")
    	  return(f)
    })
  lai.trend.chg <- do.call("rbind", ft)
    lai.sums <- tapply(lai.trend.chg$count, lai.trend.chg$metric, sum) 
  lai.trend.chg <- lai.trend.chg[grep(paste(c("-1 to -0.9", "-0.9 to -0.8", "-0.8 to -0.7", "0.7 to 0.8", 
         "0.8 to 0.9", "0.9 to 1"), collapse = "|"), lai.trend.chg$effect_size),]  
  lai.trend.chg <- lai.trend.chg[grep(paste(c("post is positive", "post is negative"), collapse = "|"), lai.trend.chg$trend),]     
  
fcov.chg <- rast("fcov_tau_change.tif")
  d <- na.omit(data.frame(c(fcov, fcov.chg)))
    ft <- lapply(c("current", "tsa"), \(i) { 
      f <- data.frame(metric = i, table(d[,"fcov_change_type"], d[,i]))
        names(f) <- c("metric",  "trend", "effect_size", "count")
    	  return(f)
    })
  fcov.trend.chg <- do.call("rbind", ft)
    fcov.sums <- tapply(fcov.trend.chg$count, fcov.trend.chg$metric, sum) 
  fcov.trend.chg <- fcov.trend.chg[grep(paste(c("-1 to -0.9", "-0.9 to -0.8", "-0.8 to -0.7", "0.7 to 0.8", 
         "0.8 to 0.9", "0.9 to 1"), collapse = "|"), fcov.trend.chg$effect_size),]  
  fcov.trend.chg <- fcov.trend.chg[grep(paste(c("post is positive", "post is negative"), collapse = "|"), fcov.trend.chg$trend),]     
