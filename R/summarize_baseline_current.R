suppressMessages(lapply(c("sf", "spatialEco", "terra", "ggplot2"), 
		         require, character.only = TRUE))

cty <- c("brazil", "bhutan", "canada", "colombia", "costa_rica", "peru")[6]

root = file.path("C:/evans/PFP/results", cty)
setwd(root)

table.dir = file.path(root, "tables", "baesline_current")
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
# Set breaks and colors
clrs <- c("darkred","aquamarine3","darkgreen")

lai.bks <- seq(0, 9.5, 0.5)
  lai.bks <- pretty_cut(lai.bks, lai.bks)[-1]
    lai.clr <- colorRampPalette(clrs)(nlevels(lai.bks))
  	  lai.labs <- levels(lai.bks)
    lai.col.table <- data.frame(value=1:length(lai.bks), col=lai.clr)				
    lai.rclass <- do.call(rbind, lapply(strsplit(as.character(lai.bks), split="to"), \(i) {
      matrix(c(as.numeric(i[1]), as.numeric(i[2])), nrow=1, ncol=2)
    }))
      lai.rclass <- cbind(lai.rclass, 1:nrow(lai.rclass))

fcov.bks <- seq(0, 1, 0.10)
  fcov.bks <- pretty_cut(fcov.bks, fcov.bks)[-1]
    fcov.clr <- colorRampPalette(clrs)(nlevels(fcov.bks))
  	  fcov.labs <- levels(fcov.bks)
    fcov.col.table <- data.frame(value=1:length(fcov.bks), col=fcov.clr)				
    fcov.rclass <- do.call(rbind, lapply(strsplit(as.character(fcov.bks), split="to"), \(i) {
      matrix(c(as.numeric(i[1]), as.numeric(i[2])), nrow=1, ncol=2)
    }))
      fcov.rclass <- cbind(fcov.rclass, 1:nrow(fcov.rclass))

#******************************************************
# Read and classify rasters
r <- list.files(root, "tif$") 
  r <- r[grep("treatment_responses", r)]

lai <- rast(r[grep("lai", r)])
  lai <- lai[[c("pre", "current")]]

  lai.class <- rast(lapply(1:nlyr(lai), \(i) {
    r <- classify(lai[[i]], rcl = lai.rclass, right=NA)
      v <- unique(r)[,1]
        levels(r) <- data.frame(value = v, class = levels(lai.bks)[v])
  	    coltab(r) <- data.frame(value = v, col = lai.clr[v])
  	return(r)
  }))
  names(lai.class) <- c("baseline", "current")
    writeRaster(lai.class, "lai_baseline_current.tif", 
                overwrite=TRUE, datatype="INT1U")
 
fcov <- rast(r[grep("fcov", r)])
  fcov <- fcov[[c("pre", "current")]]

  fcov.class <- rast(lapply(1:nlyr(fcov), \(i) {
    r <- classify(fcov[[i]], rcl = fcov.rclass, right=NA)
      v <- unique(r)[,1]
        levels(r) <- data.frame(value = v, class = levels(fcov.bks)[v])
  	    coltab(r) <- data.frame(value = v, col = fcov.clr[v])
  	return(r)
  }))
  names(fcov.class) <- c("baseline", "current")
    writeRaster(fcov.class, "fcov_baseline_current.tif", 
                overwrite=TRUE, datatype="INT1U")
 
f.lai <- freq(lai.class, usenames = TRUE)
  names(f.lai)[1:2] <- c("metric", "class")
  miss <- f.lai[FALSE,]
  f.lai$metric <- paste0("lai_", f.lai$metric) 
  f.lai$ha <- expanse(lai.class, byValue=TRUE, usenames=TRUE, unit="ha")$area
  f.lai$pct <- NA
  for(i in unique(f.lai$metric)) {
    idx <- which(f.lai$metric %in% i)
    f.lai$pct[idx] <- f.lai$count[idx] / sum(f.lai$count[idx])  
    if(length(f.lai$class[idx]) != length(lai.labs))
	  miss <- rbind(miss, data.frame(metric = i, class = 
	                lai.labs[which(!lai.labs %in% f.lai$class[idx])],
                    count = 0, ha = 0, pct = 0))			 
  }	
if(nrow(miss) > 0) f.lai <- rbind(f.lai, miss)

f.lai <- f.lai[order(f.lai$metric, f.lai$class),]
  
f.fcov <- freq(fcov.class, usenames = TRUE)
  names(f.fcov)[1:2] <- c("metric", "class")
  miss <- f.fcov[FALSE,]
  f.fcov$metric <- paste0("fcov_", f.fcov$metric) 
  f.fcov$ha <- expanse(fcov.class, byValue=TRUE, usenames=TRUE, unit="ha")$area
  f.fcov$pct <- NA
  for(i in unique(f.fcov$metric)) {
    idx <- which(f.fcov$metric %in% i)
    f.fcov$pct[idx] <- f.fcov$count[idx] / sum(f.fcov$count[idx])  
    if(length(f.fcov$class[idx]) != length(fcov.labs))
	  miss <- rbind(miss, data.frame(metric = i, class = 
	                fcov.labs[which(!fcov.labs %in% f.fcov$class[idx])],
                    count = 0, ha = 0, pct = 0))			 
  }	
if(nrow(miss) > 0) f.fcov <- rbind(f.fcov, miss)

f.fcov <- f.fcov[order(f.fcov$metric, f.fcov$class),]

es.class.freq <- rbind(f.lai, f.fcov)
  write.csv(es.class.freq, file.path(table.dir, "baseline_current_class_freq.csv"))
