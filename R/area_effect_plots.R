suppressMessages(
  lapply(c("sf", "spatialEco", "terra",  "dplyr", "ggplot2", 
           "duckdb", "DBI", "ggridges"), 
			require, character.only = TRUE))
sf::sf_use_s2(FALSE)

cty <- c("brazil", "bhutan", "canada", "colombia", "costa_rica", "peru")[5]
  root = file.path("C:/evans/PFP/results", cty)
    setwd(root)
      out.dir = file.path("C:/evans/PFP/results", cty)
      dat.dir = file.path("C:/evans/PFP", cty, "data")
    
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

lai <- rast( file.path(out.dir, "lai_effect_sizes.tif") )

m <- rast(file.path(dat.dir, "forest_500m.tif"))
  m[m < 1] <- NA
fcov <- mask(c(rast(file.path(out.dir, "fcov_effect_sizes.tif")),
             rast(file.path(out.dir, "fcov_tau.tif"), lyr=c(4,3))), m)
  names(fcov)[4:5] <- c("pre_tau", "post_tau")

  ct <- extract(fcov[[1]], interventions)
    fcov.cts <- table(ct$ID)
      fcov.pa <- interventions[-which(fcov.cts < 6),]
	    a <- units::set_units(st_area(fcov.pa), "ha") 
remove(ct, fcov.cts)

e <- extract(fcov, fcov.pa)
  for(i in 2:6) {
    e[,i] <- ifelse(e[,i] >= 0, 1, ifelse(e[,i] < 0, 0, NA))
  }
 
current.prop <- lapply(unique(e$ID), \(i) {
  p <- e[e$ID == i,][,"current"]
    p <- c(length(na.omit(p[p==0])),  
           length(na.omit(p[p==1])))
	p <- p / sum(p)
	  p[is.nan(p)] <- 0
    return(p)
})
current.prop <- do.call(rbind, current.prop)
  current.prop <- data.frame(country = cty, 
                             negative = current.prop[,1],
							 positive = current.prop[,2],
							 area = units::drop_units(a))
na.idx <- which(current.prop$negative == 0 & current.prop$positive == 0) 
  if(length(na.idx) > 0) current.prop <- current.prop[-na.idx,] 			
			
tsa.prop <- lapply(unique(e$ID), \(i) {
  p <- e[e$ID == i,][,"tsa"]
    p <- c(length(na.omit(p[p==0])),  
           length(na.omit(p[p==1])))
	p <- p / sum(p)
	  p[is.nan(p)] <- 0
    return(p)
})
tsa.prop <- do.call(rbind, tsa.prop)
  tsa.prop <- data.frame(country = cty, 
                             negative = tsa.prop[,1],
							 positive = tsa.prop[,2],
							 area = units::drop_units(a))
na.idx <- which(tsa.prop$negative == 0 & tsa.prop$positive == 0) 
  if(length(na.idx) > 0) tsa.prop <- tsa.prop[-na.idx,] 			
  
chg.prop <- lapply(unique(e$ID), \(i) {
  p <- e[e$ID == i,][,"change"]
    p <- c(length(na.omit(p[p==0])),  
           length(na.omit(p[p==1])))
	p <- p / sum(p)
	  p[is.nan(p)] <- 0
    return(p)
})
chg.prop <- do.call(rbind, chg.prop)
  chg.prop <- data.frame(country = cty, 
                             negative = chg.prop[,1],
							 positive = chg.prop[,2],
							 area = units::drop_units(a))
na.idx <- which(chg.prop$negative == 0 & chg.prop$positive == 0) 
  if(length(na.idx) > 0) chg.prop <- chg.prop[-na.idx,] 			

 
post.prop <- lapply(unique(e$ID), \(i) {
  p <- e[e$ID == i,][,"post_tau"]
    p <- c(length(na.omit(p[p==0])),  
           length(na.omit(p[p==1])))
	p <- p / sum(p)
	  p[is.nan(p)] <- 0
    return(p)
})
post.prop <- do.call(rbind, post.prop)
  post.prop <- data.frame(country = cty, 
                             negative = post.prop[,1],
							 positive = post.prop[,2],
							 area = units::drop_units(a))
na.idx <- which(post.prop$negative == 0 & post.prop$positive == 0) 
  if(length(na.idx) > 0) post.prop <- post.prop[-na.idx,] 			

par(mfrow=c(2,2))
  plot(log(post.prop$area), post.prop$positive, type="p", pch=20, col="blue",
      ylab="proportion", xlab="log(area)",
  	main="Post fcov pos/neg trend by area")
    points(log(post.prop$area), post.prop$negative, pch=20, col="red")							 
  legend("bottomleft", legend=c("positive","negative"), pch=c(20,20),
         col=c("blue", "red"), bg="white")
  plot(log(current.prop$area), current.prop$positive, type="p", pch=20, col="blue",
      ylab="proportion", xlab="log(area)",
  	main="fcov current pos/neg effect size by area")
    points(log(current.prop$area), current.prop$negative, pch=20, col="red")							 
  legend("bottomleft", legend=c("positive","negative"), pch=c(20,20),
         col=c("blue", "red"), bg="white")
  plot(log(tsa.prop$area), tsa.prop$positive, type="p", pch=20, col="blue",
      ylab="proportion", xlab="log(area)",
  	main="fcov tsa pos/neg effect size by area")
    points(log(tsa.prop$area), tsa.prop$negative, pch=20, col="red")							  
  legend("bottomleft", legend=c("positive","negative"), pch=c(20,20),
         col=c("blue", "red"), bg="white")
  plot(log(chg.prop$area), chg.prop$positive, type="p", pch=20, col="blue",
      ylab="proportion", xlab="log(area)",
  	main="fcov change pos/neg effect size by area")
    points(log(chg.prop$area), chg.prop$negative, pch=20, col="red")
  legend("bottomleft", legend=c("positive","negative"), pch=c(20,20),
         col=c("blue", "red"), bg="white")  