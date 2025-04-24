suppressMessages(lapply(c("sf", "spatialEco", "terra", "ggplot2", "geodata", 
                 "rstac","dplyr", "duckdb", "DBI", "mgcv", "ks"), 
		         require, character.only = TRUE))
sf_use_s2(FALSE)

#### which PFP to process
( cty <- c("brazil", "bhutan", "canada", "colombia", "costa_rica", "peru")[2] )

root = file.path("C:/evans/PFP", cty, "data")
setwd(root)
  dat.dir = root 
  results.dir = file.path("C:/evans/PFP/results", cty)
  mdl.dir = file.path("C:/evans/PFP", cty, "model")
  out.dir = file.path(results.dir, "climate")
    dir.create(out.dir, showWarnings = FALSE)

trends = c("tau", "slope")[2]
  switch(trends, 
    "tau" = { trend.idx = 1 },  
    "slope" = { trend.idx = 2 })

agg.by = c("intervention", "elu")[2]
vars = c("elevation", "all", "ai")[3]
scatter = c(TRUE,FALSE)[2]

dpct = 1 # deviance threshold 
thin.data = c(TRUE,FALSE)[2]
p.dat = c(0.25, 0.50, 0.75)[1] # density interval to thin
dem.res <- c(30, 90)[1] 

#*************************************************
# Process topo
#*************************************************
# if(!file.exists("elev.tif")) {
#   olm <- read_stac("http://s3.eu-central-1.wasabisys.com/stac/openlandmap/catalog.json")
#     olm$links <- links(olm, rel == "child")
#        glc_link <- links(olm, grepl("DEM", title))[[1]]
#       glc <- link_open(glc_link)
#     glc_items <- read_items(glc, progress = FALSE)
#   urls <- assets_url(glc_items, asset_names = "dsm_glo30_m_30m_s", append_gdalvsi = TRUE)
#   
#   bdy <- st_buffer(bdy, 500)
#     ref <- rast(ext(bdy), resolution = dem.res, crs=crs(bdy))
#       ref <- setValues(ref, rep(1, ncell(ref)) )
#         ref <- crop(ref, bdy)
#  
#   e <- ext(st_transform(bdy, st_crs(4326)))
#   if(dem.res == 30) {
#     elev <- rast(urls, vsi = TRUE, raw = TRUE) 
#   } else if(dem.res == 90) {
#    # add download geodata::elevation_30s() 
#   }
#     elev <- crop(elev, e)
#       elev <- mask(project(elev, ref, method = "bilinear"), bdy)
#         names(elev) <- "elevation" 
#   writeRaster(elev, "elev.tif", gdal=c("COMPRESS=LZW"), 
#               datatype="FLT4S", overwrite = TRUE)
# }			

if(vars == "elevation") {
  topo <- rast("elev.tif")
    names(topo) <- "elevation"
} else if(vars == "all"){
  elev <- rast("elev.tif")
    names(elev) <- "elevation"
  sa <- terrain(elev, v=c("slope", "aspect"), unit="degrees")
  topo <- c(elev, sa)
} else if(vars == "ai"){
  topo <- rast(file.path(dirname(root), "climate", "aridity_index.tif"))
    names(topo) <- "aridity_index"
}

#*************************************************
# Read model data
#*************************************************
lai.db = file.path(mdl.dir, "lai", "lai_250m_model_data.duckdb")
fcov.db = file.path(mdl.dir, "fcov", "fcov_500m_model_data.duckdb")
lai.results = c(file.path(results.dir, "lai_effect_sizes.tif"),
                file.path(mdl.dir, "lai/tau/lai_pre_tau_250m.tif"), 
                file.path(mdl.dir, "lai/tau/lai_post_tau_250m.tif"))
fcov.results = c(file.path(results.dir, "fcov_effect_sizes.tif"),
                 file.path(mdl.dir, "fcov/tau/fcov_pre_tau_500m.tif"), 	
                 file.path(mdl.dir, "fcov/tau/fcov_post_tau_500m.tif"))

# Read data
con <- dbConnect(duckdb::duckdb(), lai.db, read_only=TRUE)
  db.table <- dbListTables(con)
    dat <- dbGetQuery(con, paste0("SELECT * FROM ", db.table[1]))
      if(cty == "brazil") {
        dat <- rbind(dat, dbGetQuery(con, paste0("SELECT * FROM ", db.table[2])))
      }  
dbDisconnect(con, shutdown = TRUE)
  dat <- dat[-which(dat$NAME == "control"),]
  lai <- c(rast(lai.results[1]),
           rast(lai.results[2], lyr=trend.idx),
           rast(lai.results[3], lyr=trend.idx))
    names(lai)[4:5] <- c("pre", "post")

# Extract data	
dat <- data.frame(dat, elevation=extract(topo, dat[,c("X", "Y")])[,-1])
  names(dat)[grep("elevation", names(dat))] <- names(topo)
    if(vars == "all") {
      dat$hli <- 1 - cos((3.142/180)  * (dat$aspect - 30))  / 2 
      dat$scosa <- sa.trans(dat$slope, dat$aspect)
    }
	
dat <- data.frame(dat[,c("NAME", "DESIG", "IUCN", names(topo), "X", "Y")], extract(lai, dat[,c("X", "Y")])[,-1])
  if(trends == "slope"){
    dat$pre <- dat$pre / max(dat$pre, na.rm=TRUE)
    dat$post <- dat$post / max(dat$post, na.rm=TRUE)
  }

if(agg.by == "elu") {
  elu <- rast("elu.tif")
    names(elu) <- "elu"
  dat <- data.frame(dat, elu=extract(elu, dat[,c("X", "Y")])[,-1])
}

#*************************************************
# calculate areas
#*************************************************
if(agg.by == "intervention") {
  bdy <- st_read(paste0(cty, ".gpkg"), "boundary")
  if(any(cty %in% c("colombia", "peru", "brazil"))){
    pa <- st_cast(st_read(file.path(dat.dir, paste0(cty, ".gpkg")), "interventions"), "POLYGON")
      iidx <- grep("Indigenous", pa$DESIG) 
  	  if(length(iidx) > 0) pa <- pa[-iidx,]
  } else {
     pa <- st_cast(st_read(file.path(dat.dir, paste0(cty, ".gpkg")), "protected_areas"), "POLYGON")
  }
  pa.ha <- unlist(lapply(unique(pa$NAME), \(i) {
    sum(units::drop_units(units::set_units(st_area(pa[pa$NAME == i,]), "ha"))) 
  }))
  areas <- data.frame(intervention = unique(pa$NAME), area = pa.ha)
} else if(agg.by == "elu") {
  if(any(cty %in% c("colombia", "peru", "brazil"))){
    pa <- st_cast(st_read(file.path(dat.dir, paste0(cty, ".gpkg")), "interventions"), "POLYGON")
  } else {
    pa <- st_cast(st_read(file.path(dat.dir, paste0(cty, ".gpkg")), "protected_areas"), "POLYGON")
  }  
  m <- rast("forest_250m.tif")
    m[m < 1] <- NA
	  m <- mask(m, pa)
	  elu <- mask(elu, m)
  a <- expanse(elu, unit="ha", byValue=TRUE, transform=FALSE)[,-1]
    names(a)[1] <- "elu"
	  areas <- data.frame(elu=unique(a$elu), area=NA)
        for(j in unique(a$elu)) {
          areas[which(areas$elu %in% j),][2] <- sum(a[which(a$elu %in% j),]$area)
		  
        } 
}

#*************************************************
# loop for plots and table
#*************************************************
# var.combine <- unfactor(expand.grid(unique(dat$NAME), c("pre", "post"), c("elevation")))\
if(agg.by == "intervention") {
  var.combine <- unfactor(as.data.frame(table(dat$NAME)))
    names(var.combine) <- c("intervention", "freq")  
     rm.idx <- which(var.combine$freq < 20)
       if(length(rm.idx) > 0) var.combine <- var.combine[-rm.idx,]
} else if(agg.by == "elu") {
  var.combine <- unfactor(as.data.frame(table(dat$elu)))
    names(var.combine) <- c("elu", "freq")  
     rm.idx <- which(var.combine$freq < 20)
       if(length(rm.idx) > 0) var.combine <- var.combine[-rm.idx,]
}	 	 
responses <- data.frame(intervention=var.combine[,1], pre.range=NA, post.range=NA, xrange=NA, pre.slope=NA, post.slope=NA, 
                        pre.direction=NA, post.direction=NA, n=NA, pre.pValue=NA, post.pValue=NA, pre.se=NA, post.se=NA,
						welch_t=NA)
  if(agg.by == "elu") names(responses)[1] <- "elu"						

m = names(topo)[1]					
scatter.plots <- list()	  
 for(i in 1:nrow(var.combine)){
   j = as.character(var.combine[i,][1]) 
     cat("Fitting GLM between", "pre/post lai trend and elevation for", j, "(", i, "in", nrow(var.combine), ")", "\n")
        flush.console(); Sys.sleep(0.01)		
      if(agg.by == "elu") {
	    xy <- cbind(dat[dat$elu == j,][c("pre", "post", m)])
      } else if(agg.by == "intervention") {
	    xy <- cbind(dat[dat$NAME == j,][c("pre", "post", m)])
	  }
	  colnames(xy) <- c("pre", "post", m)
          na.idx <- as.data.frame(which(is.na(xy), arr.ind = TRUE))[,1]
            if(length(na.idx) > 0) xy <- xy[-na.idx,]	      
	    if(nrow(xy) > 20) {
		  ct = nrow(xy)
		  pre.mdl <- MASS::rlm(xy[,"pre"] ~ xy[,m], model=FALSE)
		  post.mdl <- MASS::rlm(xy[,"post"] ~ xy[,m], model=FALSE)
            if(scatter) {		  
              plot(xy[,m], xy[,"post"], pch=20, cex = 0.7, col="grey",   
                ylim=range(xy$pre,xy$post), xlab = "elevation", ylab = "post lai trend",  
	            main = paste0("post LAI trend vs ", m, " for ", "\n", j) )
	             abline(pre.mdl, col = "blue")
                 abline(post.mdl, col = "green")
                   legend("bottomleft", legend=c("pre slope", "post slope"), 
	            	        bg="white", lty=c(1,1), col=c("blue", "green"))
              scatter.plots[[i]] <- recordPlot()	
            }			
            responses[i,][2] <- paste0(round(min(xy[,1],na.rm=TRUE),4), " to ", round(max(xy[,1],na.rm=TRUE),4)) 
            responses[i,][3] <- paste0(round(min(xy[,2],na.rm=TRUE),4), " to ", round(max(xy[,2],na.rm=TRUE),4)) 					
            responses[i,][4] <- paste0(round(min(xy[,3],na.rm=TRUE),4), " to ", round(max(xy[,3],na.rm=TRUE),4)) 					
            responses[i,][5] <- round(coefficients(pre.mdl)[2], 6)
            responses[i,][6] <- round(coefficients(post.mdl)[2], 6)
            responses[i,][7] <- as.character(factor(sign(coefficients(pre.mdl)[2]), levels = c(-1, 0, 1), labels = c("negative", "zero", "positive"))	 )    
            responses[i,][8] <- as.character(factor(sign(coefficients(post.mdl)[2]), levels = c(-1, 0, 1), labels = c("negative", "zero", "positive"))	 )    
            responses[i,][9] <- ct	
            responses[i,][10] <- round(sfsmisc::f.robftest(pre.mdl)$p.value, 5)	
            responses[i,][11] <- round(sfsmisc::f.robftest(post.mdl)$p.value, 5)
            responses[i,][12] <- round(as.numeric(coef(summary(pre.mdl))[,2][2]),6)
            responses[i,][13] <- round(as.numeric(coef(summary(post.mdl))[,2][2]),6)
			responses[i,][14] <- (responses[i,][5] - responses[i,][6]) / sqrt(responses[i,][12]^2 + responses[i,][13]^2) 
        remove(xy,pre.mdl,post.mdl)			
      } else {
	    ct = nrow(xy)
		if(scatter) scatter.plots[[i]] <- NULL
        responses[i,][2] <- paste0(round(min(xy[,1],na.rm=TRUE),4), " to ", round(max(xy[,1],na.rm=TRUE),4)) 
        responses[i,][3] <- paste0(round(min(xy[,2],na.rm=TRUE),4), " to ", round(max(xy[,2],na.rm=TRUE),4)) 
        responses[i,][3] <- paste0(round(min(xy[,3],na.rm=TRUE),4), " to ", round(max(xy[,3],na.rm=TRUE),4)) 
        responses[i,][9] <- ct			
      }	
	}
	
	# write table
    if(agg.by == "intervention"){ 
	  responses <- merge(responses, areas, by="intervention")	
	    responses <- responses[with(responses, order(intervention)),]	
    } else if(agg.by == "elu"){
      responses <- merge(responses, areas, by="elu")		
	    responses <- responses[with(responses, order(elu)),]
    }
	responses <- data.frame(pfp=cty, responses) 		  
	  na.idx <- which(is.na(responses$pre.slope))
	    if(length(na.idx) > 0) responses <- responses[-na.idx,]
	write.csv(responses, file.path(out.dir, paste0(cty, "_lai_", vars, "_", agg.by, ".csv")), 
			            row.names = FALSE)

	# write plots
	if(scatter) {
      scatter.plots <- scatter.plots[which(!sapply(scatter.plots, is.null))]
        graph.out = file.path(out.dir, paste0(cty, "_lai_", vars, "_scatter_", agg.by, ".pdf"))
          pdf(graph.out, height=8.5, width=11)	
            for(p in 1:length(scatter.plots)) replayPlot(scatter.plots[[p]]) 
          dev.off()
    }
	
# scatter.plots <- list()	  
#   for(i in 1:nrow(var.combine)){
#       xc = as.character(var.combine[i,][3])
#       m = as.character(var.combine[i,][2])
#       j = as.character(var.combine[i,][1]) 
#       cat("Fitting GLM between", m, "and", xc, "for", j, "(", i, "in", nrow(var.combine), ")", "\n")
#         flush.console(); Sys.sleep(0.01)		
#         if(any(m %in% c("pre","post"))) {
#           lab = paste0("LAI ", m, " trend vs ", xc, " for ", j)
#         } else {
#           lab = paste0("LAI ", m, " effect size vs ", xc, " for ", j)
#         }
#       xy <- cbind(dat[dat$NAME == j,][c(m,xc)])
#         colnames(xy) <- c(m,xc)
#           na.idx <- as.data.frame(which(is.na(xy), arr.ind = TRUE))[,1]
#             if(length(na.idx) > 0) xy <- xy[-na.idx,]	  
# 		      if(thin.data & nrow(xy) > 20) {
#                 fhat <- ks::kde(xy[,1:2], h=Hpi(xy[,1:2]))
#                   p = ks::contourLevels(fhat, prob=p.dat)
#                     est <- kde(xy[,1:2], h=Hpi(xy[,1:2]), eval.points = xy[,1:2])$estimate 
# 	                  xy <- as.data.frame(xy)[which(est > p),]
#               }
# 	    if(nrow(xy) > 20) {
# 			ct = nrow(xy)
# 			  pidx <- which(xy[,m] >= 0)
# 			  nidx <- which(xy[,m] < 0)
# 			    mdl <- MASS::rlm(xy[,m] ~ xy[,xc], model=FALSE)
# 		          if(length(pidx) >= 20) 
# 			        mdl.pos <- MASS::rlm(xy[,m][pidx] ~ xy[,xc][pidx], model=FALSE)
#                   if(length(nidx) >= 20) 			 
# 			        mdl.neg <- MASS::rlm(xy[,m][nidx] ~ xy[,xc][nidx], model=FALSE)      
#             cols <- ifelse(xy[,m] >= 0, "blue", 
#                       ifelse(xy[,m] < 0, "red", NA))
# 			plot(xy[,xc], xy[,m], pch=20, col = cols, 
# 			       xlab = xc, ylab = m, main = lab)
# 				abline(mdl, col = "grey")
# 				if(exists("mdl.pos")) abline(mdl.pos, lty=3, col="blue")
# 				if(exists("mdl.neg")) abline(mdl.neg, lty=3, col="red")
#                 legend("bottomleft", legend=c("y positive", "y negative", "slope all", 
# 				       "slope pos", "slope neg"), bg="white", pch=c(20,20,NA,NA,NA), 
# 					   lty=c(NA,NA,1,3,3), col=c("blue", "red", "grey", "blue", "red"))
# 		 	scatter.plots[[i]] <- recordPlot()
# 	        responses[i,][4] <- paste0(round(min(xy[,2],na.rm=TRUE),4), " to ", round(max(xy[,2],na.rm=TRUE),4)) 
# 	        responses[i,][5] <- paste0(round(min(xy[,1],na.rm=TRUE),4), " to ", round(max(xy[,1],na.rm=TRUE),4)) 					
#             responses[i,][6] <- as.character(factor(sign(coefficients(mdl)[2]), levels = c(-1, 0, 1), labels = c("negative", "zero", "positive"))	 )    
# 		    responses[i,][7] <- coefficients(mdl)[2]
# 		    responses[i,][8] <- ct
#         remove(xy,pidx,nidx,mdl,mdl.pos,mdl.neg)			
#       } else {
# 	    ct = nrow(xy)
# 	    responses[i,][4] <- paste0(round(min(xy[,2],na.rm=TRUE),4), " to ", round(max(xy[,2],na.rm=TRUE),4)) 
# 	    responses[i,][5] <- paste0(round(min(xy[,1],na.rm=TRUE),4), " to ", round(max(xy[,1],na.rm=TRUE),4)) 
#         responses[i,][8] <- ct		
#         scatter.plots[[i]] <- NULL
#       }	
# 	}
# 	
# 	responses <- merge(responses, areas, by="intervention")	
# 	  na.idx <- which(is.na(responses$slope))
# 	    if(length(na.idx) > 0) responses <- responses[-na.idx,]
#           responses <- responses[with(responses, order(intervention, response, indicator)),]		
# 			write.csv(responses, file.path(out.dir, paste0(cty, "_lai_topo", ".csv")), row.names = FALSE)
# 	scatter.plots <- scatter.plots[which(!sapply(scatter.plots, is.null))]
# 	  graph.out = file.path(out.dir, paste0(cty, "_lai_topo_scatter.pdf"))
#         pdf(graph.out, height=8.5, width=11)	
# 	      for(p in 1:length(scatter.plots)) replayPlot(scatter.plots[[p]]) 
# 	    dev.off()

	   
# 	
# for(i in 1:nrow(var.combine)){
#    res.plots <- list()
#    nl.plots <- list()
#    scatter.plots <- list()
#      xc = as.character(var.combine[i,][3])
#      m = as.character(var.combine[i,][2])
#      j = as.character(var.combine[i,][1]) 
#      cat("Fitting GLM between", m, "and", xc, "for", j, "(", i, "in", nrow(var.combine), ")", "\n")
#        flush.console(); Sys.sleep(0.01)		
#        if(any(m %in% c("pre","post"))) {
#          lab = paste0("LAI ", m, " trend vs ", xc, " for ", j)
#        } else {
#          lab = paste0("LAI ", m, " effect size vs ", xc, " for ", j)
#        }
#      xy <- cbind(dat[dat$NAME == j,][c(m,xc,"X","Y")])
#        colnames(xy) <- c(m,xc,"X","Y")
#          na.idx <- as.data.frame(which(is.na(xy), arr.ind = TRUE))[,1]
#            if(length(na.idx) > 0) xy <- xy[-na.idx,]
#    
#	  if(nrow(xy) > 10) {
#            fhat <- ks::kde(xy[,3:4], h=Hpi(xy[,3:4]))
#              p = ks::contourLevels(fhat, prob=c(0.50))
#                est <- kde(xy[,3:4], h=Hpi(xy[,3:4]), eval.points = xy[,3:4])$estimate
#            	  if(length(na.idx) > 0) est <- insert.values(est, NA, na.idx) 
#	                xy <- as.data.frame(xy)[which(est > p),]
#
#      g <- tryCatch({ mgcv::gam(xy[,m] ~ s(xy[,xc])) }, 
#	    error = function(e) { NULL })
#		if(!is.null(g)) {	
#    	     if(round(summary(g)$dev.expl,4)*100 <  dpct) {
#               nl.plots[[i]] <- NULL
#               scatter.plots[[i]] <- NULL
#		     } else {
#              plot(g, pch=20, cex=0.1, se = TRUE, rug=FALSE, residuals=FALSE,  
#                   main=lab, xlab = xc, ylab = m, 
#    	    	      sub=paste0("R-square = ", round(summary(g)$r.sq,3),
#                    " Deviance explained = ", round(summary(g)$dev.expl,4)*100,"%")) 
#		 	nl.plots[[i]] <- recordPlot()
#              plot(g, pch=20, cex=0.1, se = TRUE, rug=FALSE, residuals=TRUE,  
#                   main=lab, xlab = xc, ylab = m, 
#    	    	      sub=paste0("R-square = ", round(summary(g)$r.sq,3),
#                    " Deviance explained = ", round(summary(g)$dev.expl,4)*100,"%")) 
#		 	res.plots[[i]] <- recordPlot()
#              mdl <- MASS::rlm(as.formula(paste0(m, "~", xc)), data=xy)			
#		 	  plot(xy[,xc], xy[,m], pch=20, cex=0.3, xlab = xc, ylab = m, main = lab)
#                 abline(mdl, col = "red")
#		 	scatter.plots[[i]] <- recordPlot()
#	        responses[i,][4] <- paste0(round(min(xy[,2],na.rm=TRUE),4), " - ", round(max(xy[,2],na.rm=TRUE),4)) 
#	        responses[i,][5] <- paste0(round(min(xy[,1],na.rm=TRUE),4), " - ", round(max(xy[,1],na.rm=TRUE),4)) 					
#            responses[i,][6] <- as.character(factor(sign(coefficients(mdl)[2]), levels = c(-1, 0, 1), labels = c("negative", "zero", "positive"))	 )    
#		    responses[i,][7] <- coefficients(mdl)[2]     
#          } 
#    	}
#      }
#      if(nrow(xy) < 10 | is.null(g)) {	
#          nl.plots[[i]] <- NULL
#          scatter.plots[[i]] <- NULL
#        }	
#	}
#	
#    write.csv(responses, file.path(out.dir, paste0(cty, "_lai_topo", ".csv")))
#	
#	nl.plots <- nl.plots[which(!sapply(nl.plots, is.null))]
#	scatter.plots <- scatter.plots[which(!sapply(scatter.plots, is.null))]
#	res.plots <- res.plots[which(!sapply(res.plots, is.null))]
#	  graph.out = file.path(out.dir, paste0(cty, "_lai_topo.pdf"))
#        pdf(graph.out, height=8.5, width=11)	
#	      for(p in 1:length(nl.plots)) replayPlot(nl.plots[[p]]) 
#	    dev.off()
#	  graph.out = file.path(out.dir, paste0(cty, "_lai_topo_scatter.pdf"))
#        pdf(graph.out, height=8.5, width=11)	
#	      for(p in 1:length(scatter.plots)) replayPlot(scatter.plots[[p]]) 
#	    dev.off()	   
#	  graph.out = file.path(out.dir, paste0(cty, "_lai_topo_residuals.pdf"))
#        pdf(graph.out, height=8.5, width=11)	
#	      for(p in 1:length(res.plots)) replayPlot(res.plots[[p]]) 
#	    dev.off()	   	  
