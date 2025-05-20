# Classifies effect size and trend rasters, returns summary tables
#  run mask.by = "interventions" first so that everything is
#  calculated then subsets (eg., IUCN) can be derived
#
# Inputs: 
# .../country/model/lai/metric/metric_250m_effect_size.tif
# .../country/model/fcov/metric/metric_500m_effect_size.tif
# .../country/model/fcov/tau/fcov_pre_tau_500m.tif
# .../country/model/fcov/tau/fcov_post_tau_500m.tif
# .../country/model/lai/tau/lai_pre_tau_250m.tif
# .../country/model/lai/tau/lai_post_tau_250m.tif
#
# Outputs: 
# Fractional Cover rasters
#   .../country/fcov_effect_sizes.tif  -  compiled fcov effect size rasters for cpct, current, tsa 
#   .../country/fcov_effect_sizes_classified.tif – Above effect size rasters classified into 21 classes (0.1 increments)
#   .../country/fcov_tau.tif – Compiled pre/post fcov trends
#   .../country/fcov_tau_classified.tif – Above trend rasters classified into 21 classes (0.1 increments)
#   .../country/fcov_tau_change.tif – classified amount of change and change direction
# 
# Leaf Area Index (LAI) rasters
#   .../country/lai_effect_sizes.tif -  compiled lai effect size rasters for cpct, current, tsa 
#   .../country/lai_effect_sizes_classified.tif – Above effect size rasters classified into 21 classes (0.1 increments)
#   .../country/lai_tau.tif– Compiled pre/post lai trends
# .../country/lai_tau_classified.tif– Above trend rasters classified into 21 classes (0.1 increments)  
# .../country/lai_tau_change.tif – classified amount of change and change direction
#   
# 
# Treatment polygons
#   .../country/ATE_polygons.gpkg – Average Treatment Effects for intervention polygons
# 
# Query tables for; all interventions, IUCN (I & II), only PA's and only PFP's
#   .../country/tables/query/effect_size_class_freq.csv – comma separated frequency, hectares and percent of classified effect size classes
#   .../country/tables/query/trend_class_freq.csv – comma separated frequency, hectares and percent of classified trend classes
#   .../country/tables/query/trend_change_freq.csv – comma separated frequency, hectares and percent of classified magnitude of change in trend
#   .../country/tables/query/trend_change_direction_freq.csv -  – comma separated frequency, hectares and percent of direction fo change in pre/post trends
# 
# The “queries” are in directories denoting their subset and the above tables exists for each subset. 
#   intervention - are all of the intervention polygons
#   pa – Are just Protected Areas excepting Bhutan, Canada and Costa Rica because the pa and intervention units are the same
#   pfp – are just the PFP’s but, only for Brazil, Colombia and Peru.  
#   Iucn – this is subset to IUCN Ia, Ib and II
#
suppressMessages(lapply(c("sf", "spatialEco", "terra", "ggplot2"), 
		         require, character.only = TRUE))
				 
#tmpFiles(current=TRUE, orphan=TRUE, old=TRUE, remove=TRUE)
#rm(list=ls(all=TRUE))
#  gc()

#### which PFP to process
( cty <- c("brazil", "bhutan", "canada", "colombia", "costa_rica", "peru")[5] )
root = file.path("C:/evans/PFP/results", cty)
setwd(root)
  out.dir = file.path("C:/evans/PFP/results", cty)
  dat.dir = file.path("C:/evans/PFP", cty, "data")
  mdl.dir = file.path("C:/evans/PFP", cty, "model")

#### processing switches
mask.by <- c("interventions", "iucn", "pfp", "pa", "country")[1] 
cal.areas = c(TRUE, FALSE)[2]
cal.es = c(TRUE, FALSE)[1]
cal.trend = c(TRUE, FALSE)[1]

#### classification paramters
n.class = c(7, 9, 21)[3]                      # number of classes
fcov.ha.scalar = 25                           # count scalar for 250m lai (average cell size in ha)
lai.ha.scalar = 6.25	                      # count scalar for 5000m fcov (average cell size in ha)
iucn.level <- c("I", "I and II")[1]
  if(cty == "peru") iucn.level = "I and II"  

lai.out <- file.path(out.dir, "lai_effect_sizes.tif")	  
fcov.out <- file.path(out.dir, "fcov_effect_sizes.tif")	  
table.dir = file.path(out.dir, "tables", mask.by)
  dir.create(table.dir, showWarnings = FALSE)

#******************************************************
# Functions
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

#******************************************************
# Set breaks and colors, note this is where we set the
# threshold around zero (and in AssignEffectSizesPoints.R)
#clrs <- c("red","yellow","cyan","blue","green")
clrs <- c("darkred","darkorange1","azure3","aquamarine3","darkgreen")
  if(n.class == 7) {  
    bks <- c(-1.00, -0.75, -0.25, -0.05, 0.05,  0.25,  0.75, 1) 
      hte.bks <- pretty_cut(bks, bks)[-1]	 
        hte.mcls <- colorRampPalette(clrs)(nlevels(hte.bks))
          hte.labs <- c("High negative", "Moderate negative", "Low negative",
		                "No effect",  "Low positive",  "Moderate positive",
						"High positive")	
  } else if(n.class == 9) {  
    bks <- c(-1.00, -0.75, -0.50, -0.25, -0.01,  0.01,  0.25,  0.50,  0.75, 1) 
      hte.bks <- pretty_cut(bks, bks)[-1]	 
        hte.mcls <- colorRampPalette(clrs)(nlevels(hte.bks))
        hte.labs <- c("High negative", "Med high negative",
        	          "Med negative", "Low negative", "No effect",
        	          "Low positive", "Med positive", "Med high positive",
        	          "High positive")	
  } else if(n.class == 21) {
    bks <- c( seq(-1,-0.005, 0.10), c(-0.01, 0.01), seq(0.10, 1, 0.1)) 
      hte.bks <- pretty_cut(bks, bks)[-1]	 
        hte.mcls <- colorRampPalette(clrs)(nlevels(hte.bks))
  	    hte.labs <- levels(hte.bks)
  }								
col.table <- data.frame(value=1:length(hte.bks), col=hte.mcls)				

# bks.clr <- read.csv("C:/evans/PFP/results/code/ColorMaps.csv")
# col.table$col %>% scales::show_col(cex_label = 0.7) 

rclass <- do.call(rbind, lapply(strsplit(as.character(hte.bks), split="to"), \(i) {
  matrix(c(as.numeric(i[1]), as.numeric(i[2])), nrow=1, ncol=2)
}))
  rclass <- cbind(rclass, 1:nrow(rclass))

#******************************************************
# Read polygon data and build PA/PFP predicates
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
    interventions <- st_read(file.path(dat.dir, paste0(cty, ".gpkg")), "interventions") 
} else {
  pa <- st_read(file.path(dat.dir, paste0(cty, ".gpkg")), "protected_areas")
    st_geometry(pa) <- "geometry"
  interventions <- pa	
}

if(iucn.level == "I") {
  iucn.idx <- which(pa$IUCN %in% "Ia" | pa$IUCN %in% "Ib")
} else if(iucn.level == "I and II") {
  iucn.idx <- which(pa$IUCN %in% "Ia" | pa$IUCN %in% "Ib" | pa$IUCN %in% "II")
}

# create IUCN aggregrated polygons (collapse Ia & Ib)
iucn <- pa[iucn.idx,]
  iucn$IUCN[which(iucn$IUCN %in% "Ia" | iucn$IUCN %in% "Ib")] <- "I"
    iucn <- iucn %>%
      dplyr::group_by(IUCN) %>% 
        dplyr::summarize(geometry = st_union(geometry))

# set mask subset polygons
switch(mask.by, 
  "country" = { pmask = bdy },  
  "pa" = { pmask = pa },  
  "pfp" = { pmask = pfp },
  "interventions" = { pmask = interventions },
  "iucn" = { pmask = pa[iucn.idx,] } )

#******************************************************
# Calculate areas for intervention and forests
if(cal.areas) {
  f250 <- rast(file.path(dat.dir, "forest_250m.tif"))
    f250[f250 == 0] <- NA
  
  imask <- rasterize(pmask, f250, background = NA)
  fmask <- f250
  cmask <- rasterize(bdy, f250, background = NA) 
    total.area <- global(cellSize(cmask, mask = TRUE, unit="ha"), sum, na.rm=TRUE)[,1]
    pa.area <- global(cellSize(imask, mask = TRUE, unit="ha"), sum, na.rm=TRUE)[,1]
    forest.cty.area <- global(mask(cellSize(fmask, mask = TRUE, unit="ha"), cmask) , sum, na.rm=TRUE)[,1] 
    forest.pa.area <- global(mask(cellSize(fmask, mask = TRUE, unit="ha"), imask), sum, na.rm=TRUE)[,1]
  
  pct.pa <- pa.area / total.area  
  pct.forest <- forest.cty.area / total.area    
  pct.pa.forest <- forest.pa.area / pa.area  
  
  areas <- data.frame(type = mask.by, country_area = total.area, 
                      forest_area = forest.cty.area, 
                      intervention_area = pa.area, 
  					forest_intervention_area = forest.pa.area,
  					intervention_pct = pct.pa, forest_pct = pct.forest,
  					forest_intervention_pct = pct.pa.forest)
  
  write.table(areas, file.path(dirname(table.dir), "areas.csv"),  
              sep = ",", row.names = FALSE, append = TRUE, 
              col.names = !file.exists(file.path(dirname(table.dir), "areas.csv")))
}
 
#*******************************************
# process effect size
#*******************************************
if(cal.es) {
#*******************************************
#*******************************************
# process effect size rasters
r <- list.files(mdl.dir, "tif$", full.names = TRUE, recursive = TRUE)
  r <- r[-grep(paste(c("tau", "gain", "loss"), collapse="|"), r)]  
    rm.idx <- grep("responses", r)
	  if(length(rm.idx) > 0) r <- r[-rm.idx] 
  
if(!file.exists(lai.out) & mask.by == "interventions"){
  lai <- rast(lapply(r[grep("250m", r)], rast, lyr = 1))
    names(lai) <- unlist(lapply(strsplit(rm.ext(basename(r[grep("250m", r)])), "_"), \(i) i[1] ))	  
      writeRaster(lai, lai.out, overwrite=TRUE, datatype="FLT4S")
} 
lai <- rast(lai.out)    
  if(any(mask.by == c("iucn", "pa", "pfp")))
    lai <- mask(lai, pmask)

if(!file.exists(fcov.out) & mask.by == "interventions"){
  fcov <- rast(lapply(r[grep("500m", r)], rast, lyr = 1))
    names(fcov) <- unlist(lapply(strsplit(rm.ext(basename(r[grep("500m", r)])), "_"), \(i) i[1] ))	
      writeRaster(fcov, fcov.out, overwrite=TRUE, datatype="FLT4S")	  
} 
fcov <- rast(fcov.out)
  if(any(mask.by == c("iucn", "pa", "pfp")))
    fcov <- mask(fcov, pmask)

#*******************************************		  
# classify lai and fcov effect size rasters
if(!file.exists(file.path(out.dir, paste0(rm.ext(basename(lai.out)), "_classified.tif")))  & mask.by == "interventions"){
  lai.class <- rast(lapply(1:nlyr(lai), \(i) {
    r <- classify(lai[[i]], rcl=rclass, right=NA)
      v <- unique(r)[,1]
        levels(r) <- data.frame(value = v, class = levels(hte.bks)[v])
  	    # coltab(r) <- data.frame(value = v, col = hte.mcls[v])
  	return(r)
  }))
  names(lai.class) <- names(lai)
    writeRaster(lai.class, file.path(out.dir, paste0(rm.ext(basename(lai.out)), "_classified.tif")), 
                overwrite=TRUE, datatype="INT1U")
}
lai.class <- rast(file.path(out.dir, paste0(rm.ext(basename(lai.out)), "_classified.tif")))
  if(any(mask.by == c("iucn", "pa", "pfp")))
    lai.class <- mask(lai.class, pmask)

if(!file.exists(file.path(out.dir, paste0(rm.ext(basename(fcov.out)), "_classified.tif"))) & mask.by == "interventions"){
  fcov.class <- rast(lapply(1:nlyr(fcov), \(i) {
    r <- classify(fcov[[i]], rcl=rclass, right=NA)
      v <- unique(r)[,1]
        levels(r) <- data.frame(value = v, class = levels(hte.bks)[v])
  	    # coltab(r) <- data.frame(value = v, col = hte.mcls[v])
  	return(r)
  }))
  names(fcov.class) <- names(fcov)
    writeRaster(fcov.class, file.path(out.dir, paste0(rm.ext(basename(fcov.out)), "_classified.tif")), 
                overwrite=TRUE, datatype="INT1U") 		
}
fcov.class <- rast(file.path(out.dir, paste0(rm.ext(basename(fcov.out)), "_classified.tif")))
  if(any(mask.by == c("iucn", "pa", "pfp")))
    fcov.class <- mask(fcov.class, pmask)

#*******************************************
# Write effect size frequency tables
f.lai <- freq(lai.class, usenames = TRUE)
  names(f.lai)[1:2] <- c("metric", "class")
  miss <- f.lai[FALSE,]
  f.lai$metric <- paste0("lai_", f.lai$metric) 
  # f.lai$ha <- f.lai$count * lai.ha.scalar
  f.lai$ha <- expanse(lai.class, byValue=TRUE, transform=FALSE, usenames=TRUE, unit="ha")$area
  f.lai$pct <- NA
  for(i in unique(f.lai$metric)) {
    idx <- which(f.lai$metric %in% i)
    f.lai$pct[idx] <- f.lai$count[idx] / sum(f.lai$count[idx])  
    if(length(f.lai$class[idx]) != length(hte.labs))
	  miss <- rbind(miss, data.frame(metric = i, class = 
	                hte.labs[which(!hte.labs %in% f.lai$class[idx])],
                    count = 0, ha = 0, pct = 0))			 
  }	
if(nrow(miss) > 0) {
  f.lai <- rbind(f.lai, miss)
  f.lai <- lapply(unique(f.lai$metric), \(f) {
    u <- f.lai[f.lai$metric == f,] 
	return(u[match(hte.labs,u$class),])
  })
  f.lai <- do.call("rbind", f.lai) 
}

f.fcov <- freq(fcov.class, usenames = TRUE)
  names(f.fcov)[1:2] <- c("metric", "class")
  miss <- f.fcov[FALSE,]
  f.fcov$metric <- paste0("fcov_", f.fcov$metric) 
  f.fcov$ha <- expanse(fcov.class, byValue=TRUE, transform=FALSE, usenames=TRUE, unit="ha")$area
  f.fcov$pct <- NA
  for(i in unique(f.fcov$metric)) {
    idx <- which(f.fcov$metric %in% i)
    f.fcov$pct[idx] <- f.fcov$count[idx] / sum(f.fcov$count[idx])  
    if(length(f.fcov$class[idx]) != length(hte.labs))
	  miss <- rbind(miss, data.frame(metric = i, class = 
	                hte.labs[which(!hte.labs %in% f.fcov$class[idx])],
                    count = 0, ha = 0, pct = 0))			 
  }	
if(nrow(miss) > 0) {
  f.fcov <- rbind(f.fcov, miss)
  f.fcov <- lapply(unique(f.fcov$metric), \(f) {
    u <- f.fcov[f.fcov$metric == f,] 
	return(u[match(hte.labs,u$class),])
  })
  f.fcov <- do.call("rbind", f.fcov) 
}

es.class.freq <- rbind(f.lai, f.fcov)
  es.class.freq <- es.class.freq[with(es.class.freq, order(metric)),]
    write.csv(es.class.freq, file.path(table.dir, "effect_size_class_freq.csv"))

} # end of effect size processing

#*******************************************
# process trend 
#*******************************************
if(cal.trend) {
#*******************************************
#*******************************************

# process effect size rasters
r <- list.files(file.path(mdl.dir, "tau"), "tif$", full.names = TRUE, recursive = TRUE)
  r <- r[grep("tau", r)]  

if(!file.exists(file.path(out.dir, "lai_tau.tif")) & mask.by == "country"){	
  lai.tau <- rast(lapply(r[grep("250m", r)], rast, lyr = 1))
    names(lai.tau) <- paste0("full_tau_", unlist(lapply(strsplit(rm.ext(basename(r[grep("250m", r)])), "_"), \(i) i[2] )))	  
      lai.tau <- c(lai.tau, lai.tau)   
	    names(lai.tau)[3:4] <- paste0("pa_tau_", unlist(lapply(strsplit(rm.ext(basename(r[grep("250m", r)])), "_"), \(i) i[2] )))	  
	      lai.tau[[3:4]] <- lai.tau[[3:4]]
      writeRaster(lai.tau, file.path(out.dir, "lai_tau.tif"), overwrite=TRUE, datatype="FLT4S")  
}
lai.tau <- rast(file.path(out.dir, "lai_tau.tif"), lyr=c(1,2))
  names(lai.tau)<- gsub("full_", "", names(lai.tau))
if(any(mask.by == c("interventions", "iucn", "pa", "pfp")))
  lai.tau <- mask(lai.tau, pmask)
  
if(!file.exists(file.path(out.dir, "fcov_tau.tif")) & mask.by == "country"){
  fcov.tau <- rast(lapply(r[grep("500m", r)], rast, lyr = 1))
    names(fcov.tau) <- paste0("full_tau_", unlist(lapply(strsplit(rm.ext(basename(r[grep("500m", r)])), "_"), \(i) i[2] )))	
      fcov.tau <- c(fcov.tau, fcov.tau)   
	    names(fcov.tau)[3:4] <- paste0("pa_tau_", unlist(lapply(strsplit(rm.ext(basename(r[grep("500m", r)])), "_"), \(i) i[2] )))	  
	      fcov.tau[[3:4]] <- mask(fcov.tau[[3:4]], pmask)  
      writeRaster(fcov.tau, file.path(out.dir, "fcov_tau.tif"), overwrite=TRUE, datatype="FLT4S") 
}
fcov.tau <- rast(file.path(out.dir, "fcov_tau.tif"), lyr=c(1,2))
  names(fcov.tau)<- gsub("full_", "", names(fcov.tau))    
if(any(mask.by == c("interventions", "iucn", "pa", "pfp")))
  fcov.tau <- mask(fcov.tau, pmask)

#*******************************************		  
# classify trend rasters
if(!file.exists(file.path(out.dir, "lai_tau_classified.tif")) & mask.by == "country"){
  lai.tau.class <- rast(lapply(1:nlyr(lai.tau), \(i) {
    r <- classify(lai.tau[[i]], rcl=rclass, right=NA)
      v <- unique(r)[,1]
        levels(r) <- data.frame(value = v, class = levels(hte.bks)[v])
  	    # coltab(r) <- data.frame(value = v, col = hte.mcls[v])
  	return(r)
  }))
  names(lai.tau.class) <- names(lai.tau)
    writeRaster(lai.tau.class, file.path(out.dir, "lai_tau_classified.tif"), 
                overwrite=TRUE, datatype="FLT4S")		
}
lai.tau.class <- rast(file.path(out.dir, "lai_tau_classified.tif"), lyr=c(1,2))
  names(lai.tau.class)<- gsub("full_", "", names(lai.tau.class))
if(any(mask.by == c("interventions", "iucn", "pa", "pfp")))
  lai.tau.class <- mask(lai.tau.class, pmask)
 
if(!file.exists(file.path(out.dir, "fcov_tau_classified.tif")) & mask.by == "country"){
  fcov.tau.class <- rast(lapply(1:nlyr(fcov.tau), \(i) {
    r <- classify(fcov.tau[[i]], rcl=rclass, right=NA)
      v <- unique(r)[,1]
        levels(r) <- data.frame(value = v, class = levels(hte.bks)[v])
  	    # coltab(r) <- data.frame(value = v, col = hte.mcls[v])
  	return(r)
  }))
  names(fcov.tau.class) <- names(fcov.tau)
    writeRaster(fcov.tau.class, file.path(out.dir, "fcov_tau_classified.tif"), 
                overwrite=TRUE, datatype="FLT4S")
}
fcov.tau.class <- rast(file.path(out.dir, "fcov_tau_classified.tif"), lyr=c(1,2))
  names(fcov.tau.class)<- gsub("full_", "", names(fcov.tau.class))
if(any(mask.by == c("interventions", "iucn", "pa", "pfp")))
  fcov.tau.class <- mask(fcov.tau.class, pmask)

#*******************************************
# Write trend frequency tables
f.lai <- freq(lai.tau.class, usenames = TRUE)
  names(f.lai)[1:2] <- c("period", "class")
  f.lai$period <- paste0("lai_", f.lai$period) 
  f.lai$ha <- expanse(lai.tau.class, byValue=TRUE, transform=FALSE, usenames=TRUE, unit="ha")$area
  #f.lai$ha <- f.lai$count * lai.ha.scalar
  f.lai$pct <- NA
  miss <- f.lai[FALSE,]
    for(i in unique(f.lai$period)) {
      idx <- which(f.lai$period %in% i)
      f.lai$pct[idx] <- f.lai$count[idx] / sum(f.lai$count[idx])
      if(length(f.lai$class[idx]) != length(hte.labs)){
	    miss <- rbind(miss, data.frame(period = i, class = 
	                  hte.labs[which(!hte.labs %in% f.lai$class[idx])],
                      count = 0, ha = 0, pct = 0))
	  }			  
    }
if(nrow(miss) > 0) f.lai <- rbind(f.lai, miss)
	
f.fcov <- freq(fcov.tau.class, usenames = TRUE)
  names(f.fcov)[1:2] <- c("period", "class")
  f.fcov$period <- paste0("fcov_", f.fcov$period) 
  f.fcov$ha <- expanse(fcov.tau.class, byValue=TRUE, transform=FALSE, usenames=TRUE, unit="ha")$area
  #f.fcov$ha <-  f.fcov$count * fcov.ha.scalar
  f.fcov$pct <- NA
  miss <- f.fcov[FALSE,]
    for(i in unique(f.fcov$period)) {
      idx <- which(f.fcov$period %in% i)
      f.fcov$pct[idx] <- f.fcov$count[idx] / sum(f.fcov$count[idx])
      if(length(f.fcov$class[idx]) != length(hte.labs)) {
	    miss <- rbind(miss, data.frame(period = i, class = 
	                  hte.labs[which(!hte.labs %in% f.fcov$class[idx])],
                      count = 0, ha = 0, pct = 0))
      }					  
    }
if(nrow(miss) > 0) f.fcov <- rbind(f.fcov, miss)
	
tau.class.freq <- rbind(f.lai, f.fcov)
  write.csv(tau.class.freq, file.path(table.dir, "trend_class_freq.csv"))

#******************************************************
# classify direction of change in trends
#******************************************************

#**** Classify direction of change
sign.class <- function(x) {
  s <- sign(x)
  ifelse(s[1] == -1 & s[2] == -1, 1,
    ifelse(s[1] == -1 & s[2] == 1, 2,
      ifelse(s[1] == 1 & s[2] == -1, 3,
	    ifelse(s[1] == 1 & s[2] == 1, 4,
      ifelse(s[1] == 0 & s[2] == 0, 5,
    ifelse(s[1] == 0 & s[2] == -1, 6,
  ifelse(s[1] == -1 & s[2] == 0, 7, NA)))))))
}

chg.classes <- data.frame(values=c(1:7),
  class= c("pre is negative - post is negative", 
           "pre is negative - post is positive",
           "pre is positive - post is negative",
		   "pre is positive - post is positive",
		   "pre is no change - post is no change",		   
		   "pre is no change - post is negative",
		   "pre is negative - post is no change"))

if(!file.exists(file.path(out.dir, "lai_tau_change.tif")) & mask.by == "country"){
  lai.chg <- app(c(lai.tau[[2]], lai.tau[[1]]), \(x) { sign.class(x) } )
    v <- unique(lai.chg)[,1]
      levels(lai.chg) <- chg.classes[which(chg.classes$values %in% v),] 
        names(lai.chg) <- "lai_change_type"
    writeRaster(lai.chg, file.path(out.dir, "lai_tau_change.tif"), overwrite=TRUE)
} 
lai.chg <- rast(file.path(out.dir, "lai_tau_change.tif"))
if(any(mask.by == c("interventions", "iucn", "pa", "pfp")))
  lai.chg <- mask(lai.chg, pmask)

if(!file.exists(file.path(out.dir, "fcov_tau_change.tif")) & mask.by == "country"){
  fcov.chg <- app(c(fcov.tau[[2]], fcov.tau[[1]]), \(x) { sign.class(x) } )
    v <- unique(fcov.chg)[,1]
      levels(fcov.chg) <- chg.classes[which(chg.classes$values %in% v),] 
       names(fcov.chg) <- "fcov_change_type"
    writeRaster(fcov.chg, file.path(out.dir, "fcov_tau_change.tif"), overwrite=TRUE) 
}
fcov.chg <- rast(file.path(out.dir, "fcov_tau_change.tif"))
if(any(mask.by == c("interventions", "iucn", "pa", "pfp")))
  fcov.chg <- mask(fcov.chg, pmask)

#*******************************************
# Write trend frequency tables

# Trend change direction frequencies
lai.dir <- freq(lai.chg[[1]], usenames = TRUE)
  names(lai.dir)[1:2] <- c("metric", "class")
    lai.dir$ha <- expanse(lai.chg[[1]], byValue=TRUE, transform=FALSE, usenames=TRUE, unit="ha")$area
    lai.dir$pct <- lai.dir$count / sum(lai.dir$count)   
	miss <- lai.dir[FALSE,]
      if(length(unique(lai.dir$class)) != length(chg.classes$class)) {	
        miss <- rbind(miss, data.frame(metric = "lai_change_type", class = 
  	                  chg.classes$class[which(!chg.classes$class %in% lai.dir$class)],
                      count = 0, ha = 0, pct = 0))
      }
if(nrow(miss) > 0) lai.dir <- rbind(lai.dir, miss)

fcov.dir <- freq(fcov.chg[[1]], usenames = TRUE)
  names(fcov.dir)[1:2] <- c("metric", "class")
    fcov.dir$ha <- expanse(fcov.chg[[1]], byValue=TRUE, transform=FALSE, usenames=TRUE, unit="ha")$area
    fcov.dir$pct <- fcov.dir$count / sum(fcov.dir$count)   
	miss <- fcov.dir[FALSE,]
      if(length(unique(fcov.dir$class)) != length(chg.classes$class)) {	
        miss <- rbind(miss, data.frame(metric = "fcov_change_type", class = 
  	                  chg.classes$class[which(!chg.classes$class %in% fcov.dir$class)],
                      count = 0, ha = 0, pct = 0))
      }
if(nrow(miss) > 0) fcov.dir <- rbind(fcov.dir, miss)
	 
tau.dir <- rbind(lai.dir, fcov.dir) 
write.csv(tau.dir, file.path(table.dir, "slope_change_direction_freq.csv"))

#*******************************************
# Trend direction by effect size classes

if(mask.by != "country"){
  d <- na.omit(data.frame(c(lai.class, lai.chg)))
    ft <- lapply(c("change", "current", "tsa"), \(i) { 
      f <- data.frame(metric = paste0("lai_", i), table(d[,"lai_change_type"], d[,i]))
        names(f) <- c("metric",  "trend", "effect_size", "count")
    	return(f)
    })
  lai.f <- do.call("rbind", ft)
   
  d <- na.omit(data.frame(c(fcov.class, fcov.chg)))
    ft <- lapply(c("change", "current", "tsa"), \(i) { 
      f <- data.frame(metric = paste0("fcov_", i), table(d[,"fcov_change_type"], d[,i]))
        names(f) <- c("metric",  "trend", "effect_size", "count")
    	return(f)
    })
  fcov.f <- do.call("rbind", ft)
  
  f <- rbind(lai.f, fcov.f)
  
  write.csv(f, file.path(table.dir, "effect_size_trend_frequencies.csv"), 
            row.names = FALSE)
  }

} # end of trend processing

#*******************************************
#*******************************************
## assign ATE to intervention polygons
# if(!file.exists(file.path(out.dir, "ATE_polygons.gpkg"))){
#   if(any(cty %in% c("colombia", "peru", "brazil"))){ 	  
#     pa <- st_read(file.path(dat.dir, paste0(cty, ".gpkg")), "protected_areas")
#       st_geometry(pa) <- "geometry" 
#         if(cty == "brazil") pa <- pa[-which(pa$DESIG %in% "Indigenous Area"),]
#     	  pa$type <- "protected area"	
#       elai <- extract(lai, pa)
#         dlai <- data.frame(lai_cpct = as.numeric(tapply(elai$cpct, elai$ID, median, na.rm=TRUE)),
#                            lai_current =  as.numeric(tapply(elai$current, elai$ID, median, na.rm=TRUE)), 
#                            lai_tsa = as.numeric(tapply(elai$tsa, elai$ID, median, na.rm=TRUE) ))
#       efcov <- extract(fcov, pa)
#         dfcov <- data.frame(fcov_cpct = as.numeric(tapply(efcov$cpct, efcov$ID, median, na.rm=TRUE)),
#                             fcov_current =  as.numeric(tapply(efcov$current, efcov$ID, median, na.rm=TRUE)), 
#                             fcov_tsa = as.numeric(tapply(efcov$tsa, efcov$ID, median, na.rm=TRUE) ))
#             d <- cbind(dlai, dfcov)
#           pa <- cbind(pa, d)
#       st_write(pa, file.path(out.dir, "ATE_polygons.gpkg"), "protected_areas", append=FALSE)
#     
# 	pfp <- st_read(file.path(dat.dir, paste0(cty, ".gpkg")), "PFP")
#       st_geometry(pfp) <- "geometry"
#     	pfp$IUCN <- "Not Reported"
#   	      pfp$type <- "PFP"
#       elai <- extract(lai, pfp)
#         dlai <- data.frame(lai_cpct = as.numeric(tapply(elai$cpct, elai$ID, median, na.rm=TRUE)),
#                            lai_current =  as.numeric(tapply(elai$current, elai$ID, median, na.rm=TRUE)), 
#                            lai_tsa = as.numeric(tapply(elai$tsa, elai$ID, median, na.rm=TRUE) ))
#       efcov <- extract(fcov, pfp)
#         dfcov <- data.frame(fcov_cpct = as.numeric(tapply(efcov$cpct, efcov$ID, median, na.rm=TRUE)),
#                             fcov_current =  as.numeric(tapply(efcov$current, efcov$ID, median, na.rm=TRUE)), 
#                             fcov_tsa = as.numeric(tapply(efcov$tsa, efcov$ID, median, na.rm=TRUE) ))
#         d <- cbind(dlai, dfcov)
#           pfp <- cbind(pfp, d)
#       st_write(pfp, file.path(out.dir, "ATE_polygons.gpkg"), "PFP", append=FALSE)
#   } else {
#     pa <- st_read(file.path(dat.dir, paste0(cty, ".gpkg")), "protected_areas")
#       elai <- extract(lai, pa)
#         dlai <- data.frame(lai_cpct = as.numeric(tapply(elai$cpct, elai$ID, median, na.rm=TRUE)),
#                            lai_current =  as.numeric(tapply(elai$current, elai$ID, median, na.rm=TRUE)), 
#                            lai_tsa = as.numeric(tapply(elai$tsa, elai$ID, median, na.rm=TRUE) ))
#       efcov <- extract(fcov, pa)
#         dfcov <- data.frame(fcov_cpct = as.numeric(tapply(efcov$cpct, efcov$ID, median, na.rm=TRUE)),
#                             fcov_current =  as.numeric(tapply(efcov$current, efcov$ID, median, na.rm=TRUE)), 
#                             fcov_tsa = as.numeric(tapply(efcov$tsa, efcov$ID, median, na.rm=TRUE) ))
#         d <- cbind(dlai, dfcov)
#           pa <- cbind(pa, d)
#         # na.idx <- unique(which(is.na(pa), arr.ind = TRUE)[,"row"])
#         # if(length(na.idx) > 0) pa <- pa[-na.idx,]		
#         # plot(pa[c("fcov_cpct", "lai_cpct")], border=NA)
#       st_write(pa, file.path(out.dir, "ATE_polygons.gpkg"), "protected_areas", append=FALSE)
#   }
# }
# 
# 
# v <- levels(lai.class$cpct)[[1]]$value 
# coltab(lai.class$cpct) <- data.frame(value = v, col = hte.mcls[v])
# 
# v <- levels(lai.class$current)[[1]]$value 
# coltab(lai.class$current) <- data.frame(value = v, col = hte.mcls[v])
# 
# v <- levels(lai.class$tsa)[[1]]$value 
# coltab(lai.class$tsa) <- data.frame(value = v, col = hte.mcls[v])
 