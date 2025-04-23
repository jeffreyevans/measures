suppressMessages(lapply(c("sf", "spatialEco", "terra", "ggplot2"), 
		         require, character.only = TRUE))

idx = 1
( cty <- c("bhutan", "costa_rica", "canada", "colombia", "peru", "brazil")[idx] )

root = file.path("C:/evans/PFP/results", cty)
setwd(root)
  out.dir = file.path("C:/evans/PFP/results", cty, "tables")
  dat.dir = file.path("C:/evans/PFP", cty, "data")
  tau.dir = file.path("C:/evans/PFP", cty, "model", "tau")

#******************************************************
# read intervention polygons
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
    pmask <- st_read(file.path(dat.dir, paste0(cty, ".gpkg")), "intervention_boundaries") 
} else {
  pa <- st_read(file.path(dat.dir, paste0(cty, ".gpkg")), "protected_areas")
    st_geometry(pa) <- "geometry"
  pmask <- pa	
}

#******************************************************
# Set breaks and colors
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
bks <- pretty_cut(bks, bks)[-1]
mcls <- colorRampPalette(clrs)(nlevels(bks))
rclass <- do.call(rbind, lapply(strsplit(as.character(bks), split="to"), \(i) {
  matrix(c(as.numeric(i[1]), as.numeric(i[2])), nrow=1, ncol=2)
}))
  rclass <- cbind(rclass, 1:nrow(rclass))

#******************************************************
# Read forest mask and create cell area rasters
f <- list.files(dat.dir, "tif$", full.names=TRUE)
  f <- f[grep("forest", f)]
f250 <- rast(grep("250", f, value=TRUE))
  a250 <- cellSize(f250, mask=TRUE, unit="ha")
f500 <- rast(grep("500", f, value=TRUE))
  a500 <- cellSize(f250, mask=TRUE, unit="ha")

# read LAI and fCOV for entire contry
r <- list.files(tau.dir, "tif$", full.names=TRUE)
lai.tau <- rast(r[grep("lai", r)])
  names(lai.tau)[c(1,8)] <- c("post", "pre")
    lai.tau <- lai.tau[[c(8,1)]] 
fcov.tau <- rast(r[grep("fcov", r)])
  names(fcov.tau)[c(1,8)] <- c("post", "pre")
    fcov.tau <- fcov.tau[[c(8,1)]] 

#******************************************************
# Classify by tau breaks
lai.tau.class <- rast(lapply(1:nlyr(lai.tau), \(i) {
  r <- classify(lai.tau[[i]], rcl=rclass, right=NA)
    v <- unique(r)[,1]
      levels(r) <- data.frame(value = v, class = levels(bks)[v])
	    coltab(r) <- data.frame(value = v, col = mcls[v])
	return(r)
}))
names(lai.tau.class) <- names(lai.tau)
  
fcov.tau.class <- rast(lapply(1:nlyr(fcov.tau), \(i) {
   r <- classify(fcov.tau[[i]], rcl=rclass, right=NA)
     v <- unique(r)[,1]
       levels(r) <- data.frame(value = v, class = levels(bks)[v])
 	     coltab(r) <- data.frame(value = v, col = mcls[v])
 	return(r)
 }))
 names(fcov.tau.class) <- names(fcov.tau)


diff.lai$ha <- expanse(lai.chg[[1]], byValue=TRUE, usenames=TRUE, unit="ha")$area