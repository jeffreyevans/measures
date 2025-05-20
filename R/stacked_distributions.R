suppressMessages(
  lapply(c("sf", "spatialEco", "terra",  "dplyr", "ggplot2", 
           "duckdb", "DBI", "ggridges"), 
			require, character.only = TRUE))

cty <- c("brazil", "bhutan", "canada", "colombia", "costa_rica", "peru")

theme_set(theme_minimal())
clrs <- c("darkred","darkorange1","azure3","aquamarine3","darkgreen")
  bks <- c( seq(-1,-0.005, 0.10), c(-0.01, 0.01), seq(0.10, 1, 0.1)) 
    bks <- cut(bks, bks)[-1]	 
      mcls <- colorRampPalette(clrs)(nlevels(bks))

agg.by = c("ecosystem", "climate", "iucn")[3]
type = "effect_sizes"
trends = c("tau", "slope")[1]
  switch(trends, 
    "tau" = { trend.idx = 1 },  
    "slope" = { trend.idx = 2 })

for(j in cty) {
 cat("Calculating", agg.by, "effect size/change/trend distributions for", j, "\n")
   flush.console(); Sys.sleep(0.01)	
  setwd(file.path("C:/evans/PFP", j, "model")) 
  	lai.db = file.path(getwd(), "lai", "lai_250m_model_data.duckdb")
    fcov.db = file.path(getwd(), "fcov", "fcov_500m_model_data.duckdb")

    graph.out = file.path("C:/evans/PFP/results", j, paste0(j, "_", agg.by, "_distributions", ".pdf"))  
	
    lai.results = c(file.path("C:/evans/PFP/results", j, "lai_effect_sizes.tif"),
                    file.path("C:/evans/PFP", j, "model/lai/tau/lai_pre_tau_250m.tif"), 
                    file.path("C:/evans/PFP", j, "model/lai/tau/lai_post_tau_250m.tif"))
    fcov.results = c(file.path("C:/evans/PFP/results", j, "fcov_effect_sizes.tif"),
                     file.path("C:/evans/PFP", j, "model/fcov/tau/fcov_pre_tau_500m.tif"), 	
                     file.path("C:/evans/PFP", j, "model/fcov/tau/fcov_post_tau_500m.tif"))

    elu <- rast(file.path("C:/evans/PFP", j, "data", "elu.tif"))
	  names(elu) <- "ecosystem"
    clim <- rast(file.path("C:/evans/PFP", j, "data", "elu_temperature_moisture.tif"))
      names(clim) <- "bioclimatic"
  
	if(agg.by == "iucn") {
      pa <- st_read( file.path("C:/evans/PFP", j, "data", paste0(j,".gpkg")), "protected_areas")
	    pa <- pa[,-which(names(pa) %in% c("NAME", "DESIG"))]
          pa[grep(paste(c("Ia", "Ib"), collapse="|"), pa$IUCN),]$IUCN <- "Ia & Ib"   		
    }  
  
  pdf(graph.out, height=8.5, width=11)	  

    #****************************************************    
    #****************************************************
    # LAI 
    #****************************************************
    #****************************************************	
	con <- dbConnect(duckdb::duckdb(), lai.db, read_only=TRUE)
      db.table <- dbListTables(con)
      dat <- dbGetQuery(con, paste0("SELECT * FROM ", db.table[1]))
	if(j == "brazil") {
	  dat <- rbind(dat, dbGetQuery(con, paste0("SELECT * FROM ", db.table[2])))
	}  
    dbDisconnect(con, shutdown = TRUE)
	  dat <- dat[-which(dat$NAME == "control"),]
	    lai <- c(rast(lai.results[1]),
	             rast(lai.results[2], lyr=trend.idx),
	             rast(lai.results[3], lyr=trend.idx))
          names(lai)[4:5] <- c("pre", "post")
	  
	dat <- data.frame(dat, extract(c(elu, clim), dat[,c("X", "Y")])[,-1])
	  dat <- data.frame(dat[,c("NAME", "DESIG", "IUCN", "ecosystem", "bioclimatic", "X", "Y")], 
	                    extract(lai, dat[,c("X", "Y")])[,-1])
        if(trends == "slope"){
		  neg.idx <- which(dat$pre < 0)
		  pos.idx <- which(dat$pre <= 0)
		  dat$pre[neg.idx] <- dat$pre[neg.idx] / abs(min(dat$pre[neg.idx], na.rm=TRUE))
		  dat$pre[pos.idx] <- dat$pre[pos.idx] / max(dat$pre[pos.idx], na.rm=TRUE)
		  neg.idx <- which(dat$post < 0)
		  pos.idx <- which(dat$post <= 0)	      
		  dat$post[neg.idx] <- dat$post[neg.idx] / abs(min(dat$post[neg.idx], na.rm=TRUE))
		  dat$post[pos.idx] <- dat$post[pos.idx] / max(dat$post[pos.idx], na.rm=TRUE)
	    }
    #### CLIMATE ####		
	if(agg.by == "climate") {
      dat <- dat[-which(dat$bioclimatic == "converted"),] 
	    na.idx <- which(is.na(dat$bioclimatic))
	      if(length(na.idx) > 1) dat <- dat[-na.idx,] 
	  elu.freq <- table(dat$bioclimatic)
	    elu.pct <- ( elu.freq / sum(elu.freq) ) * 100
          elu.drop <- names(elu.pct)[which(elu.pct < 1)]
            if(length(elu.drop) > 0) {
		      elu.pct <- elu.pct[-which(names(elu.pct) %in% elu.drop)]
              dat <- dat[-which(dat$bioclimatic %in% elu.drop),]
            }
		  elu.pct <- data.frame(bioclimatic = names(elu.pct), percent = round(as.numeric(elu.pct),2))
		    elu.pct$change <- max(dat$change, na.rm=TRUE)
			elu.pct$current <- max(dat$current, na.rm=TRUE)
		    elu.pct$tsa <- max(dat$tsa, na.rm=TRUE)			
		    elu.pct$pre <- max(dat$pre, na.rm=TRUE)
		    elu.pct$post <- max(dat$post, na.rm=TRUE)
				
      # plot by bioclimatic temperature/moisture
      curr.plot <- ggplot(dat, aes(x = current, y = bioclimatic, group = bioclimatic, fill = after_stat(x))) +
        geom_density_ridges_gradient(scale=1, calc_ecdf = TRUE) +
          scale_fill_gradientn(colours = mcls, name = "current") +
		  geom_vline(xintercept=0, linetype="dashed", color = "black") +
            labs(title = 'LAI Current Effect Size by bioclimatic gradient') +	
			  geom_text(data = elu.pct, aes(y=bioclimatic, label = paste0(percent,"%")),
                        position=position_nudge(y = 0.5), size=3.5) 	
        print(curr.plot)	  
      tsa.plot <- ggplot(dat, aes(x = tsa, y = bioclimatic, group = bioclimatic, fill = after_stat(x))) +
        geom_density_ridges_gradient(scale=1, calc_ecdf = TRUE) +
          scale_fill_gradientn(colours = mcls, name = "tsa") +
		  geom_vline(xintercept=0, linetype="dashed", color = "black") +
            labs(title = 'LAI TSA Effect Size by bioclimatic gradient') +	
			  geom_text(data = elu.pct, aes(y=bioclimatic, label = paste0(percent,"%")),
                        position=position_nudge(y = 0.5), size=3.5) 	
        print(tsa.plot) 
      chg.plot <- ggplot(dat, aes(x = change, y = bioclimatic, group = bioclimatic, fill = after_stat(x))) +
        geom_density_ridges_gradient(scale=1, calc_ecdf = TRUE) +
          scale_fill_gradientn(colours = mcls, name = "change") +
		  geom_vline(xintercept=0, linetype="dashed", color = "black") +
            labs(title = 'LAI Absolute Change by bioclimatic gradient') +	
			  geom_text(data = elu.pct, aes(y=bioclimatic, label = paste0(percent,"%")),
                        position=position_nudge(y = 0.5), size=3.5) 	
        print(chg.plot)		
      pre.trend.plot <- ggplot(dat, aes(x = pre, y = bioclimatic, group = bioclimatic, fill = after_stat(x))) +
        geom_density_ridges_gradient(scale=1, calc_ecdf = TRUE) +
          scale_fill_gradientn(colours = mcls, name = "trend") +
		  geom_vline(xintercept=0, linetype="dashed", color = "black") +
		    scale_x_continuous(limits = range(c(dat$pre,dat$post),na.rm=TRUE)) +
              labs(title = 'LAI pre-intervention trends by bioclimatic gradient') +	
			    geom_text(data = elu.pct, aes(y=bioclimatic, label = paste0(percent,"%")),
                          position=position_nudge(y = 0.5), size=3.5)	
        print(pre.trend.plot)	
      post.trend.plot <- ggplot(dat, aes(x = post, y = bioclimatic, group = bioclimatic, fill = after_stat(x))) +
        geom_density_ridges_gradient(scale=1, calc_ecdf = TRUE) +
          scale_fill_gradientn(colours = mcls, name = "trend") +
		  geom_vline(xintercept=0, linetype="dashed", color = "black") +
		    scale_x_continuous(limits = range(c(dat$pre,dat$post),na.rm=TRUE)) +
              labs(title = 'LAI post-intervention trends by bioclimatic gradient') +	
			    geom_text(data = elu.pct, aes(y=bioclimatic, label = paste0(percent,"%")),
                          position=position_nudge(y = 0.5), size=3.5) 	
        print(post.trend.plot)	

    #### ECOSYSTEM ####
    } else if(agg.by == "ecosystem") {
     dat <- dat[-which(dat$ecosystem == "converted"),] 
	    na.idx <- which(is.na(dat$ecosystem))
	      if(length(na.idx) > 1) dat <- dat[-na.idx,] 
	  
	  elu.freq <- table(dat$ecosystem)
	    elu.pct <- ( elu.freq / sum(elu.freq) ) * 100
          elu.drop <- names(elu.pct)[which(elu.pct < 1)]
            if(length(elu.drop) > 0) {
		      elu.pct <- elu.pct[-which(names(elu.pct) %in% elu.drop)]
              dat <- dat[-which(dat$ecosystem %in% elu.drop),]
            }
		  elu.pct <- data.frame(ecosystem = names(elu.pct), percent = round(as.numeric(elu.pct),2))
		    elu.pct$current <- max(dat$current, na.rm=TRUE)
		    elu.pct$tsa <- max(dat$tsa, na.rm=TRUE)
			elu.pct$change <- max(dat$change, na.rm=TRUE)
		    elu.pct$pre <- max(dat$pre, na.rm=TRUE)
		    elu.pct$post <- max(dat$post, na.rm=TRUE)
			
      # plot by ecosystem
      curr.plot <- ggplot(dat, aes(x = current, y = ecosystem, group = ecosystem, fill = after_stat(x))) +
        geom_density_ridges_gradient(scale=1, calc_ecdf = TRUE) +
          scale_fill_gradientn(colours = mcls, name = "current") +
		  geom_vline(xintercept=0, linetype="dashed", color = "black") +
            labs(title = 'LAI Current Effect Size by ecosystem type')+ 	
			  geom_text(data = elu.pct, aes(y=ecosystem, label = paste0(percent,"%")),
                        position=position_nudge(y = 0.5), size=3.5)
        print(curr.plot)	  
      tsa.plot <- ggplot(dat, aes(x = tsa, y = ecosystem, group = ecosystem, fill = after_stat(x))) +
        geom_density_ridges_gradient(scale=1, calc_ecdf = TRUE) +
          scale_fill_gradientn(colours = mcls, name = "tsa") +
		  geom_vline(xintercept=0, linetype="dashed", color = "black") +
            labs(title = 'LAI TSA Effect Size by ecosystem type') + 	
			  geom_text(data = elu.pct, aes(y = ecosystem, label = paste0(percent,"%")),
                        position=position_nudge(y = 0.5), size=3.5)
       print(tsa.plot) 
      chg.plot <- ggplot(dat, aes(x = change, y = ecosystem, group = ecosystem, fill = after_stat(x))) +
        geom_density_ridges_gradient(scale=1, calc_ecdf = TRUE) +
          scale_fill_gradientn(colours = mcls, name = "change") +
		  geom_vline(xintercept=0, linetype="dashed", color = "black") +
            labs(title = 'LAI Absolute Change by ecosystem type')+ 	
			  geom_text(data = elu.pct, aes(y = ecosystem, label = paste0(percent,"%")),
                        position=position_nudge(y = 0.5), size=3.5)
        print(chg.plot)		
      pre.trend.plot <- ggplot(dat, aes(x = pre, y = ecosystem, group = ecosystem, fill = after_stat(x))) +
        geom_density_ridges_gradient(scale=1, calc_ecdf = TRUE) +
          scale_fill_gradientn(colours = mcls, name = "trend") +
		  geom_vline(xintercept=0, linetype="dashed", color = "black") +
		    scale_x_continuous(limits = range(c(dat$pre,dat$post),na.rm=TRUE)) +
              labs(title = 'LAI pre-intervention trends by ecosystem type')+ 	
			    geom_text(data = elu.pct, aes(y=ecosystem, label = paste0(percent,"%")),
                          position=position_nudge(y = 0.5), size=3.5)
        print(pre.trend.plot)	
      post.trend.plot <- ggplot(dat, aes(x = post, y = ecosystem, group = ecosystem, fill = after_stat(x))) +
        geom_density_ridges_gradient(scale=1, calc_ecdf = TRUE) +
          scale_fill_gradientn(colours = mcls, name = "trend") +
		  geom_vline(xintercept=0, linetype="dashed", color = "black") +
		    scale_x_continuous(limits = range(c(dat$pre,dat$post),na.rm=TRUE)) +
              labs(title = 'LAI post-intervention trends by ecosystem type') +	
			    geom_text(data = elu.pct, aes(y=ecosystem, label = paste0(percent,"%")),
                          position=position_nudge(y = 0.5), size=3.5)
        print(post.trend.plot)
 
    #### IUCN ####		
    } else if(agg.by == "iucn") {
	  dat$ID <- 1:nrow(dat)
	    dat.df <- dat
	      dat <- st_as_sf(dat[,-which(names(dat) %in% c("NAME", "DESIG", "IUCN"))], coords = c("X", "Y"), crs = st_crs(elu), agr = "constant")
	        dat <- st_drop_geometry(st_as_sf(terra::intersect(vect(dat), vect(pa))))
	  if(any(j %in% c("brazil", "colombia", "peru"))) {     
	    dat.df <- dat.df[which(!dat.df$ID %in% dat$ID),]
		  dat.df$IUCN <- "PFP"
		    dat <- rbind(dat[,-which(names(dat) %in% c("YEAR"))], 
			  dat.df[,-which(names(dat.df) %in% c("NAME", "DESIG", "X", "Y"))])
	  }
	  na.idx <- which(is.na(dat$IUCN))
	    if(length(na.idx) > 1) dat <- dat[-na.idx,]  
	  iucn.freq <- table(dat$IUCN)
	    iucn.pct <- (iucn.freq / sum(iucn.freq) ) * 100
          iucn.drop <- names(iucn.pct)[which(iucn.pct < 1)]
            if(length(iucn.drop) > 0) {
		      iucn.pct <- iucn.pct[-which(names(iucn.pct) %in% iucn.drop)]
              dat <- dat[-which(dat$IUCN %in% iucn.drop),]
            }
		  iucn.pct <- data.frame(IUCN = names(iucn.pct), percent = round(as.numeric(iucn.pct),2))
		    iucn.pct$current <- max(dat$current, na.rm=TRUE)
		    iucn.pct$tsa <- max(dat$tsa, na.rm=TRUE)
			iucn.pct$change <- max(dat$change, na.rm=TRUE)
		    iucn.pct$pre <- max(dat$pre, na.rm=TRUE)
		    iucn.pct$post <- max(dat$post, na.rm=TRUE)
	  
	  dat$add <- "all"
	  
      # plot by IUCN
      curr.plot <- ggplot(dat, aes(x = current, y = IUCN, group = IUCN, fill = after_stat(x))) +
        geom_density_ridges_gradient(scale=1, calc_ecdf = TRUE) +
          scale_fill_gradientn(colours = mcls, name = "current") +
		  geom_vline(xintercept=0, linetype="dashed", color = "black") +
            labs(title = 'LAI Current Effect Size by IUCN')+ 	
			  geom_text(data = iucn.pct, aes(y=IUCN, label = paste0(percent,"%")),
                        position=position_nudge(y = 0.5), size=3.5)
        print(curr.plot)	  
      tsa.plot <- ggplot(dat, aes(x = tsa, y = IUCN, group = IUCN, fill = after_stat(x))) +
        geom_density_ridges_gradient(scale=1, calc_ecdf = TRUE) +
          scale_fill_gradientn(colours = mcls, name = "tsa") +
		  geom_vline(xintercept=0, linetype="dashed", color = "black") +
            labs(title = 'LAI TSA Effect Size by IUCN type') + 	
			  geom_text(data = iucn.pct, aes(y = IUCN, label = paste0(percent,"%")),
                        position=position_nudge(y = 0.5), size=3.5)
       print(tsa.plot) 
      chg.plot <- ggplot(dat, aes(x = change, y = IUCN, group = IUCN, fill = after_stat(x))) +
        geom_density_ridges_gradient(scale=1, calc_ecdf = TRUE) +
          scale_fill_gradientn(colours = mcls, name = "change") +
		    geom_vline(xintercept=0, linetype="dashed", color = "black") +
              labs(title = 'LAI Absolute Change by IUCN type')+ 	
			    geom_text(data = iucn.pct, aes(y = IUCN, label = paste0(percent,"%")),
                          position=position_nudge(y = 0.5), size=3.5)
        print(chg.plot)	
      pre.trend.plot <- ggplot(dat, aes(x = pre, y = IUCN, group = IUCN, fill = after_stat(x))) +
        geom_density_ridges_gradient(scale=1, calc_ecdf = TRUE) +
          scale_fill_gradientn(colours = mcls, name = "trend") +
		  geom_vline(xintercept=0, linetype="dashed", color = "black") +
		    scale_x_continuous(limits = range(c(dat$pre,dat$post),na.rm=TRUE)) +
              labs(title = 'LAI pre-intervention trends by IUCN type')+ 	
			    geom_text(data = iucn.pct, aes(y=IUCN, label = paste0(percent,"%")),
                          position=position_nudge(y = 0.5), size=3.5)
        print(pre.trend.plot)	
      post.trend.plot <- ggplot(dat, aes(x = post, y = IUCN, group = IUCN, fill = after_stat(x))) +
        geom_density_ridges_gradient(scale=1, calc_ecdf = TRUE) +
          scale_fill_gradientn(colours = mcls, name = "trend") +
		  geom_vline(xintercept=0, linetype="dashed", color = "black") +
		    scale_x_continuous(limits = range(c(dat$pre,dat$post),na.rm=TRUE)) +
              labs(title = 'LAI post-intervention trends by IUCN type') +	
			    geom_text(data = iucn.pct, aes(y=IUCN, label = paste0(percent,"%")),
                          position=position_nudge(y = 0.5), size=3.5)
        print(post.trend.plot)	
    }	
	
    #****************************************************	
    #****************************************************
    # fCOV 
    #****************************************************
    #****************************************************	
    con <- dbConnect(duckdb::duckdb(), fcov.db, read_only=TRUE)
      db.table <- dbListTables(con)
      dat <- dbGetQuery(con, paste0("SELECT * FROM ", db.table))
    dbDisconnect(con, shutdown = TRUE)
	dat <- dat[-which(dat$NAME == "control"),]
	  fcov <- c(rast(fcov.results[1]),
	           rast(fcov.results[2], lyr=trend.idx),
	           rast(fcov.results[3], lyr=trend.idx))
          names(fcov)[4:5] <- c("pre", "post")

	dat <- data.frame(dat, extract(c(elu, clim), dat[,c("X", "Y")])[,-1])
	  dat <- data.frame(dat[,c("NAME", "DESIG", "IUCN", "ecosystem", "bioclimatic", "X", "Y")], 
	                    extract(fcov, dat[,c("X", "Y")])[,-1])
        if(trends == "slope"){
		  neg.idx <- which(dat$pre < 0)
		  pos.idx <- which(dat$pre <= 0)
		  dat$pre[neg.idx] <- dat$pre[neg.idx] / abs(min(dat$pre[neg.idx], na.rm=TRUE))
		  dat$pre[pos.idx] <- dat$pre[pos.idx] / max(dat$pre[pos.idx], na.rm=TRUE)
		  neg.idx <- which(dat$post < 0)
		  pos.idx <- which(dat$post <= 0)	      
		  dat$post[neg.idx] <- dat$post[neg.idx] / abs(min(dat$post[neg.idx], na.rm=TRUE))
		  dat$post[pos.idx] <- dat$post[pos.idx] / max(dat$post[pos.idx], na.rm=TRUE)
	    }

    #### CLIMATE ####
	if(agg.by == "climate") {    
      dat <- dat[-which(dat$bioclimatic == "converted"),] 
	    na.idx <- which(is.na(dat$bioclimatic))
	      if(length(na.idx) > 1) dat <- dat[-na.idx,] 
	  
	  elu.freq <- table(dat$bioclimatic)
	    elu.pct <- ( elu.freq / sum(elu.freq) ) * 100
          elu.drop <- names(elu.pct)[which(elu.pct < 1)]
            if(length(elu.drop) > 0) {
		      elu.pct <- elu.pct[-which(names(elu.pct) %in% elu.drop)]
              dat <- dat[-which(dat$bioclimatic %in% elu.drop),]
            }
		  elu.pct <- data.frame(bioclimatic = names(elu.pct), percent = round(as.numeric(elu.pct),2))
		    elu.pct$current <- max(dat$current, na.rm=TRUE)
		    elu.pct$tsa <- max(dat$tsa, na.rm=TRUE)	
		    elu.pct$change <- max(dat$change, na.rm=TRUE)			
		    elu.pct$pre <- max(dat$pre, na.rm=TRUE)
		    elu.pct$post <- max(dat$post, na.rm=TRUE)
				 
    # plot by temperature/moisture
    curr.plot <- ggplot(dat, aes(x = current, y = bioclimatic, group = bioclimatic, fill = after_stat(x))) +
      geom_density_ridges_gradient(scale=1, calc_ecdf = TRUE) +
        scale_fill_gradientn(colours = mcls, name = "current") +
		geom_vline(xintercept=0, linetype="dashed", color = "black") +
          labs(title = 'fCOV Current Effect Size by bioclimatic gradient') +	
			  geom_text(data = elu.pct, aes(y=bioclimatic, label = paste0(percent,"%")),
                        position=position_nudge(y = 0.5), size=3.5) 	 	
      print(curr.plot)	
    tsa.plot <- ggplot(dat, aes(x = tsa, y = bioclimatic, group = bioclimatic, fill = after_stat(x))) +
      geom_density_ridges_gradient(scale=1, calc_ecdf = TRUE) +
        scale_fill_gradientn(colours = mcls, name = "tsa") +
		geom_vline(xintercept=0, linetype="dashed", color = "black") +
          labs(title = 'fCOV TSA Effect Size by bioclimatic gradient') +	
			  geom_text(data = elu.pct, aes(y=bioclimatic, label = paste0(percent,"%")),
                        position=position_nudge(y = 0.5), size=3.5) 	 	
      print(tsa.plot)
    chg.plot <- ggplot(dat, aes(x = change, y = bioclimatic, group = bioclimatic, fill = after_stat(x))) +
      geom_density_ridges_gradient(scale=1, calc_ecdf = TRUE) +
        scale_fill_gradientn(colours = mcls, name = "change") +
		geom_vline(xintercept=0, linetype="dashed", color = "black") +
          labs(title = 'fCOV Change Effect Size by bioclimatic gradient') +	
			  geom_text(data = elu.pct, aes(y=bioclimatic, label = paste0(percent,"%")),
                        position=position_nudge(y = 0.5), size=3.5) 	 	
      print(chg.plot)	  
    pre.trend.plot <- ggplot(dat, aes(x = pre, y = bioclimatic, group = bioclimatic, fill = after_stat(x))) +
      geom_density_ridges_gradient(scale=1, calc_ecdf = TRUE) +
        scale_fill_gradientn(colours = mcls, name = "trend") +
		geom_vline(xintercept=0, linetype="dashed", color = "black") +
		  scale_x_continuous(limits = range(c(dat$pre,dat$post),na.rm=TRUE)) +
            labs(title = 'fCOV pre-intervention trends by bioclimatic gradient') +	
			  geom_text(data = elu.pct, aes(y=bioclimatic, label = paste0(percent,"%")),
                        position=position_nudge(y = 0.5), size=3.5) 	 	
      print(pre.trend.plot)	
    post.trend.plot <- ggplot(dat, aes(x = post, y = bioclimatic, group = bioclimatic, fill = after_stat(x))) +
      geom_density_ridges_gradient(scale=1, calc_ecdf = TRUE) +
        scale_fill_gradientn(colours = mcls, name = "trend") +
		geom_vline(xintercept=0, linetype="dashed", color = "black") +
		  scale_x_continuous(limits = range(c(dat$pre,dat$post),na.rm=TRUE)) +
            labs(title = 'fCOV post-intervention trends by bioclimatic gradient') +	
			  geom_text(data = elu.pct, aes(y=bioclimatic, label = paste0(percent,"%")),
                        position=position_nudge(y = 0.5), size=3.5) 		
      print(post.trend.plot)

    #### ECOSYSTEM ####	  
    } else if(agg.by == "ecosystem") {
     dat <- dat[-which(dat$ecosystem == "converted"),] 
	    na.idx <- which(is.na(dat$ecosystem))
	      if(length(na.idx) > 1) dat <- dat[-na.idx,] 
	  
	  elu.freq <- table(dat$ecosystem)
	    elu.pct <- ( elu.freq / sum(elu.freq) ) * 100
          elu.drop <- names(elu.pct)[which(elu.pct < 1)]
            if(length(elu.drop) > 0) {
		      elu.pct <- elu.pct[-which(names(elu.pct) %in% elu.drop)]
              dat <- dat[-which(dat$ecosystem %in% elu.drop),]
            }
		  elu.pct <- data.frame(ecosystem = names(elu.pct), percent = round(as.numeric(elu.pct),2))
		    elu.pct$current <- max(dat$current, na.rm=TRUE)
		    elu.pct$tsa <- max(dat$tsa, na.rm=TRUE)			
		    elu.pct$change <- max(dat$change, na.rm=TRUE)
		    elu.pct$pre <- max(dat$pre, na.rm=TRUE)
		    elu.pct$post <- max(dat$post, na.rm=TRUE)
			
    # plot by ecosystem	
    curr.plot <- ggplot(dat, aes(x = current, y = ecosystem, group = ecosystem, fill = after_stat(x))) +
      geom_density_ridges_gradient(scale=1, calc_ecdf = TRUE) +
        scale_fill_gradientn(colours = mcls, name = "current") +
		geom_vline(xintercept=0, linetype="dashed", color = "black") +
          labs(title = 'fCOV Current Effect Size by ecosystem type') +	
		    geom_text(data = elu.pct, aes(y=ecosystem, label = paste0(percent,"%")),
                      position=position_nudge(y = 0.5), size=3.5)	
      print(curr.plot)	
    tsa.plot <- ggplot(dat, aes(x = tsa, y = ecosystem, group = ecosystem, fill = after_stat(x))) +
      geom_density_ridges_gradient(scale=1, calc_ecdf = TRUE) +
        scale_fill_gradientn(colours = mcls, name = "tsa") +
		geom_vline(xintercept=0, linetype="dashed", color = "black") +
		#scale_fill_viridis_c(name = "current", direction = -1) +
          labs(title = 'fCOV TSA Effect Size by ecosystem type') +	
			    geom_text(data = elu.pct, aes(y=ecosystem, label = paste0(percent,"%")),
                          position=position_nudge(y = 0.5), size=3.5) 	
      print(tsa.plot)
    chg.plot <- ggplot(dat, aes(x = change, y = ecosystem, group = ecosystem, fill = after_stat(x))) +
      geom_density_ridges_gradient(scale=1, calc_ecdf = TRUE) +
        scale_fill_gradientn(colours = mcls, name = "change") +
		geom_vline(xintercept=0, linetype="dashed", color = "black") +
          labs(title = 'fCOV Change Effect Size by ecosystem type') +	
		    geom_text(data = elu.pct, aes(y=ecosystem, label = paste0(percent,"%")),
                      position=position_nudge(y = 0.5), size=3.5)	
      print(chg.plot)	  
    pre.trend.plot <- ggplot(dat, aes(x = pre, y = ecosystem, group = ecosystem, fill = after_stat(x))) +
      geom_density_ridges_gradient(scale=1, calc_ecdf = TRUE) +
        scale_fill_gradientn(colours = mcls, name = "trend") +
		geom_vline(xintercept=0, linetype="dashed", color = "black") +
		  scale_x_continuous(limits = range(c(dat$pre,dat$post),na.rm=TRUE)) +
            labs(title = 'fCOV pre-intervention trends by ecosystem type') +	
			  geom_text(data = elu.pct, aes(y=ecosystem, label = paste0(percent,"%")),
                        position=position_nudge(y = 0.5), size=3.5) 	
      print(pre.trend.plot)	
    post.trend.plot <- ggplot(dat, aes(x = post, y = ecosystem, group = ecosystem, fill = after_stat(x))) +
      geom_density_ridges_gradient(scale=1, calc_ecdf = TRUE) +
        scale_fill_gradientn(colours = mcls, name = "trend") +
		geom_vline(xintercept=0, linetype="dashed", color = "black") +
		  scale_x_continuous(limits = range(c(dat$pre,dat$post),na.rm=TRUE)) +
            labs(title = 'fCOV post-intervention trends by ecosystem type') +	
			  geom_text(data = elu.pct, aes(y=ecosystem, label = paste0(percent,"%")),
                        position=position_nudge(y = 0.5), size=3.5) 	
      print(post.trend.plot)	

    #### IUCN ####		
    } else if(agg.by == "iucn") {
	  dat$ID <- 1:nrow(dat)
	    dat.df <- dat
	      dat <- st_as_sf(dat[,-which(names(dat) %in% c("NAME", "DESIG", "IUCN"))], coords = c("X", "Y"), crs = st_crs(elu), agr = "constant")
	        dat <- st_drop_geometry(st_as_sf(terra::intersect(vect(dat), vect(pa))))
	  if(any(j %in% c("brazil", "colombia", "peru"))) {     
	    dat.df <- dat.df[which(!dat.df$ID %in% dat$ID),]
		  dat.df$IUCN <- "PFP"
		    dat <- rbind(dat[,-which(names(dat) %in% c("YEAR"))], 
			  dat.df[,-which(names(dat.df) %in% c("NAME", "DESIG", "X", "Y"))])
	  }
	  na.idx <- which(is.na(dat$IUCN))
	    if(length(na.idx) > 1) dat <- dat[-na.idx,]  
	  iucn.freq <- table(dat$IUCN)
	    iucn.pct <- (iucn.freq / sum(iucn.freq) ) * 100
          iucn.drop <- names(iucn.pct)[which(iucn.pct < 1)]
            if(length(iucn.drop) > 0) {
		      iucn.pct <- iucn.pct[-which(names(iucn.pct) %in% iucn.drop)]
              dat <- dat[-which(dat$IUCN %in% iucn.drop),]
            }
		  iucn.pct <- data.frame(IUCN = names(iucn.pct), percent = round(as.numeric(iucn.pct),2))
		    iucn.pct$current <- max(dat$current, na.rm=TRUE)
		    iucn.pct$tsa <- max(dat$tsa, na.rm=TRUE)
			iucn.pct$change <- max(dat$change, na.rm=TRUE)
		    iucn.pct$pre <- max(dat$pre, na.rm=TRUE)
		    iucn.pct$post <- max(dat$post, na.rm=TRUE)
			
      # plot by IUCN
      curr.plot <- ggplot(dat, aes(x = current, y = IUCN, group = IUCN, fill = after_stat(x))) +
        geom_density_ridges_gradient(scale=1, calc_ecdf = TRUE) +
          scale_fill_gradientn(colours = mcls, name = "current") +
		  geom_vline(xintercept=0, linetype="dashed", color = "black") +
            labs(title = 'FCOV Current Effect Size by IUCN')+ 	
			  geom_text(data = iucn.pct, aes(y=IUCN, label = paste0(percent,"%")),
                        position=position_nudge(y = 0.5), size=3.5)
        print(curr.plot)	  
      tsa.plot <- ggplot(dat, aes(x = tsa, y = IUCN, group = IUCN, fill = after_stat(x))) +
        geom_density_ridges_gradient(scale=1, calc_ecdf = TRUE) +
          scale_fill_gradientn(colours = mcls, name = "tsa") +
		  geom_vline(xintercept=0, linetype="dashed", color = "black") +
            labs(title = 'FCOV TSA Effect Size by IUCN type') + 	
			  geom_text(data = iucn.pct, aes(y = IUCN, label = paste0(percent,"%")),
                        position=position_nudge(y = 0.5), size=3.5)
       print(tsa.plot) 
      chg.plot <- ggplot(dat, aes(x = change, y = IUCN, group = IUCN, fill = after_stat(x))) +
        geom_density_ridges_gradient(scale=1, calc_ecdf = TRUE) +
          scale_fill_gradientn(colours = mcls, name = "change") +
		  geom_vline(xintercept=0, linetype="dashed", color = "black") +
            labs(title = 'FCOV Absolute Change by IUCN type')+ 	
			  geom_text(data = iucn.pct, aes(y = IUCN, label = paste0(percent,"%")),
                        position=position_nudge(y = 0.5), size=3.5)
        print(chg.plot)		
      pre.trend.plot <- ggplot(dat, aes(x = pre, y = IUCN, group = IUCN, fill = after_stat(x))) +
        geom_density_ridges_gradient(scale=1, calc_ecdf = TRUE) +
          scale_fill_gradientn(colours = mcls, name = "trend") +
		  geom_vline(xintercept=0, linetype="dashed", color = "black") +
		    scale_x_continuous(limits = range(c(dat$pre,dat$post),na.rm=TRUE)) +
              labs(title = 'FCOV pre-intervention trends by IUCN type')+ 	
			    geom_text(data = iucn.pct, aes(y=IUCN, label = paste0(percent,"%")),
                          position=position_nudge(y = 0.5), size=3.5)
        print(pre.trend.plot)	
      post.trend.plot <- ggplot(dat, aes(x = post, y = IUCN, group = IUCN, fill = after_stat(x))) +
        geom_density_ridges_gradient(scale=1, calc_ecdf = TRUE) +
          scale_fill_gradientn(colours = mcls, name = "trend") +
		  geom_vline(xintercept=0, linetype="dashed", color = "black") +
		    scale_x_continuous(limits = range(c(dat$pre,dat$post),na.rm=TRUE)) +
              labs(title = 'FCOV post-intervention trends by IUCN type') +	
			    geom_text(data = iucn.pct, aes(y=IUCN, label = paste0(percent,"%")),
                          position=position_nudge(y = 0.5), size=3.5)
        print(post.trend.plot)	
    }	  
  dev.off()	
}
