suppressMessages(
  lapply(c("sf", "spatialEco", "terra",  "dplyr", 
           "duckdb", "DBI"), 
			require, character.only = TRUE))

cty <- c("brazil", "bhutan", "canada", "colombia", "costa_rica", "peru")[c(4,6)]

#*************************************************
#*************************************************
# leaf area index
#*************************************************
#*************************************************
pdf(file.path("C:/evans/PFP/results", "LAI_yearly_ecosystem_trends.pdf"), height=8.5, width=11)
  for(j in cty) {
    cat("Calculating trends for", j, "\n")
    setwd(file.path("C:/evans/PFP", j, "data"))
      rn <- list.files(getwd(), "tif$") 
	  lai <- rast(grep("LAI_250m_trend", rn, value = TRUE))
        dates <- as.Date(names(lai))
          years <- format(dates, "%Y")
      f <- rast("forest_250m.tif")
        f[f < 1] <- NA	  
	  if(file.exists("elu_250m.tif")) {
	    elu <- rast("elu_250m.tif")
	  } else {
	    elu <- rast("elu_30m.tif")
		  elu <- resample(elu, f, method = "near")
	  }
	elu <- mask(elu, f)  
      elu.dat <- levels(elu)[[1]]
        levels(elu) <- NULL
		  felu <- freq(elu)
		    drop.idx <- which(felu$count < 100)
          if(length(drop.idx) > 0) felu <- felu[-drop.idx,]
      trends <- invisible(lapply(felu$value, \(e) {
	     cat("Processing ELU", e, "\n")
		   flush.console(); Sys.sleep(0.01)	
	     elu.sub <- elu
		  elu.sub[elu.sub != e] <- NA
	       trends <- vector()
	         for(i in unique(years)) {
		       cat("  calculating moments for", i, "in", max(years), "\n")
		         flush.console(); Sys.sleep(0.01)	
               r <- mask(lai[[which(years == i)]], elu.sub)
		        trends <- append(trends, mean(global(r, mean, na.rm=TRUE)[,1]))
             }
	      return(trends)
	    }))
    plot(unique(years), trends[[1]], type="n", 
	     ylim = range(unlist(lapply(trends, range))),  
	     xlab = "year", ylab = "Median LAI by Ecosystem Type",  
         main = paste0(j, " LAI  yearly trends"))
      for(i in 1:length(trends)){
	    lines(unique(years), trends[[i]])
          points(unique(years), trends[[i]], pch=20, cex = 0.75)
	  }
	abline(lm(apply(do.call(rbind, trends), MARGIN=2, mean) ~ as.numeric(unique(years))), col="red")	      
}
dev.off()


#*************************************************
# simple single trend
pdf(file.path("C:/evans/PFP/results", "LAI_yearly_trends_peru.pdf"), height=8.5, width=11)
  for(j in cty) {
    cat("Calculating trends for", j, "\n")
    setwd(file.path("C:/evans/PFP", j, "data"))
      rn <- list.files(getwd(), "tif$") 
	  lai <- rast(grep("LAI_250m_trend", rn, value = TRUE))
        dates <- as.Date(names(lai))
          years <- format(dates, "%Y")
      f <- rast("forest_250m.tif")
        f[f < 1] <- NA		  
    trends <- vector()
      for(i in unique(years)) {
	    cat("  calculating moments for", i, "in", max(years), "\n")
		  flush.console(); Sys.sleep(0.01)	
        r <- mask(lai[[which(years == i)]], f)
		trends <- append(trends, mean(global(r, mean, na.rm=TRUE)[,1]))
      }
    plot(unique(years), trends, type="o",  
	     xlab = "year", ylab = "Median LAI",  
         main = paste0(j, " LAI  yearly trends"))
      abline(lm(trends ~ as.numeric(unique(years))), col="red")
    remove(trends, lai)	
  gc()	
  }		 
dev.off()

#*************************************************
# simple trend, treatment, control, outside
pdf(file.path("C:/evans/PFP/results", "LAI_yearly_trends2.pdf"), height=8.5, width=11)
  for(j in cty) {
    cat("Calculating trends for", j, "\n")
    setwd(file.path("C:/evans/PFP", j, "data"))
	  f <- rast("forest_250m.tif")
	    f[f < 1] <- NA
	      fcells <- as.data.frame(f, cells = TRUE, na.rm = TRUE)$cell 
      rn <- list.files(getwd(), "tif$") 
	    lai <- rast(grep("lai_250m_trend", rn, value = TRUE))
          dates <- as.Date(names(lai))
            years <- format(dates, "%Y")
	r <- as.data.frame(lai, cells=TRUE, na.rm = TRUE)
	  r <- r[which(r$cell %in% fcells),]
        cells <- r[,1]
	      r <- r[,-1]
    lai.con <- file.path("C:/evans/PFP", j, "model", "lai", "lai_250m_model_data.duckdb")
      con <- dbConnect(duckdb::duckdb(), lai.con, read_only=TRUE)
       db.table <- dbListTables(con)[1]
       dat <- dbGetQuery(con, paste0("SELECT * FROM ", db.table))     
    dbDisconnect(con, shutdown = TRUE)
	ctl.idx <- dat[which(dat$NAME == "control"),]$cell
    int.idx <- dat[which(dat$NAME != "control"),]$cell
	remove(dat, fcells)
	ctl.trends <- vector()  
    int.trends <- vector()
    out.trends <- vector()
      for(i in unique(years)) {
	    cat("  calculating moments for", i, "in", max(years), "\n")
		  flush.console(); Sys.sleep(0.01)	
	    ctl.trends <- append(ctl.trends, median(as.matrix(r[which(cells %in% ctl.idx),][which(years %in% i)]), na.rm = TRUE))  
		int.trends <- append(int.trends, median(as.matrix(r[which(cells %in% int.idx),][which(years %in% i)]), na.rm = TRUE))
		out.trends <- append(out.trends, median(as.matrix(r[cells[which(!cells %in% c(ctl.idx,int.idx))],][which(years %in% i)]), na.rm = TRUE))
	  }
    plot(unique(years),int.trends, type="n", 
	     ylim=c(range(c(ctl.trends, int.trends, out.trends))),  
	     xlab = "year", ylab = "Median lai",  
         main = paste0(j, " lai  yearly trends"))
      lines(unique(years), ctl.trends, lty = 1, col = "red")
        points(unique(years), ctl.trends, pch=20, cex = 0.75, col = "red")
      lines(unique(years), int.trends, lty = 1, col = "green")
        points(unique(years), int.trends, pch=20, cex = 0.75, col = "green") 
      lines(unique(years), out.trends, lty = 3, col = "black")
        points(unique(years), out.trends, pch=20, cex = 0.75, col = "black")
	abline(lm(apply(rbind(ctl.trends, int.trends, out.trends), MARGIN=2, mean) ~ as.numeric(unique(years))), 
	       col="black", lty=6, lwd=0.5)		
    legend("bottomright", legend=c("Control", "Treatment", "Outside", "average trend"), lty=c(1,1,3, 6), 
           col=c("red", "green", "black", "black"))  
    remove(lai, ctl.idx, int.idx, ctl.trends, int.trends, out.trends, lai.con)	
  gc()	
  }		 
dev.off()


#*************************************************
#*************************************************
# fractional cover
#*************************************************
#*************************************************
pdf(file.path("C:/evans/PFP/results", "fCOV_yearly_ecosystem_trends.pdf"), height=8.5, width=11)
  for(j in cty) {
    cat("Calculating trends for", j, "\n")
    setwd(file.path("C:/evans/PFP", j, "data"))
      rn <- list.files(getwd(), "tif$") 
	  fcov <- rast(grep("fcov_500m_trend", rn, value = TRUE))
        dates <- as.Date(names(fcov))
          years <- format(dates, "%Y")
      f <- rast("forest_500m.tif")
        f[f < 1] <- NA	  
	  if(file.exists("elu_500m.tif")) {
	    elu <- rast("elu_500m.tif")
	  } else {
	    elu <- rast("elu_30m.tif")
		  elu <- resample(elu, f, method = "near")
	  }
	elu <- mask(elu, f)  
      elu.dat <- levels(elu)[[1]]
        levels(elu) <- NULL
		  felu <- freq(elu)
		    drop.idx <- which(felu$count < 100)
          if(length(drop.idx) > 0) felu <- felu[-drop.idx,]
      trends <- invisible(lapply(felu$value, \(e) {
	     cat("Processing ELU", e, "\n")
		   flush.console(); Sys.sleep(0.01)	
	     elu.sub <- elu
		  elu.sub[elu.sub != e] <- NA
	       trends <- vector()
	         for(i in unique(years)) {
		       cat("  calculating moments for", i, "in", max(years), "\n")
		         flush.console(); Sys.sleep(0.01)	
               r <- mask(fcov[[which(years == i)]], elu.sub)
		        trends <- append(trends, mean(global(r, mean, na.rm=TRUE)[,1]))
             }
	      return(trends)
	    }))
    plot(unique(years), trends[[1]], type="n", 
	     ylim = range(unlist(lapply(trends, range))),  
	     xlab = "year", ylab = "Median fCOV by Ecosystem Type",  
         main = paste0(j, " fCOV  yearly trends"))
      for(i in 1:length(trends)){
	    lines(unique(years), trends[[i]])
          points(unique(years), trends[[i]], pch=20, cex = 0.75)
	  }
	abline(lm(apply(do.call(rbind, trends), MARGIN=2, mean) ~ as.numeric(unique(years))), col="red")	      
}
dev.off()


#*************************************************
# simple single trend
pdf(file.path("C:/evans/PFP/results", "FCOV_yearly_trends_colombia.pdf"), height=8.5, width=11)
  for(j in cty) {
    cat("Calculating trends for", j, "\n")
    setwd(file.path("C:/evans/PFP", j, "data"))
      rn <- list.files(getwd(), "tif$") 
	  fcov <- rast(grep("fcov_500m_trend", rn, value = TRUE))
        dates <- as.Date(names(fcov))
          years <- format(dates, "%Y")
      f <- rast("forest_500m.tif")
        f[f < 1] <- NA		  
    trends <- vector()
      for(i in unique(years)) {
	    cat("  calculating moments for", i, "in", max(years), "\n")
		  flush.console(); Sys.sleep(0.01)	
        r <- mask(fcov[[which(years == i)]], f)
		trends <- append(trends, mean(global(r, mean, na.rm=TRUE)[,1]))
      }
    plot(unique(years), trends, type="o",  
	     xlab = "year", ylab = "Median fCOV",  
         main = paste0(j, " fCOV  yearly trends"))
      abline(lm(trends ~ as.numeric(unique(years))), col="red")
    remove(trends, fcov)	
  gc()	
  }		 
dev.off()




#*************************************************
#*************************************************
# simple trend, treatment, control, outside
pdf(file.path("C:/evans/PFP/results", "FCOV_yearly_trends.pdf"), height=8.5, width=11)
  for(j in cty) {
    cat("Calculating trends for", j, "\n")
    setwd(file.path("C:/evans/PFP", j, "data"))
	  f <- rast("forest_500m.tif")
	    f[f < 1] <- NA
	      fcells <- as.data.frame(f, cells = TRUE, na.rm = TRUE)$cell 
      rn <- list.files(getwd(), "tif$") 
	    fcov <- rast(grep("fcov_500m_trend", rn, value = TRUE))
		  #names(fcov) <- names(rast(grep("lai_250m_trend", rn, value = TRUE)))[-1005]
          dates <- as.Date(names(fcov))
            years <- format(dates, "%Y")
    fcov.con <- file.path("C:/evans/PFP", j, "model", "fcov", "fcov_500m_model_data.duckdb")
      con <- dbConnect(duckdb::duckdb(), fcov.con, read_only=TRUE)
       db.table <- dbListTables(con)[1]
       dat <- dbGetQuery(con, paste0("SELECT * FROM ", db.table))
    dbDisconnect(con, shutdown = TRUE)
	  ctl.idx <- dat[which(dat$NAME == "control"),]$cell
      int.idx <- dat[which(dat$NAME != "control"),]$cell
	  remove(dat)
	ctl.trends <- vector()  
    int.trends <- vector()
    out.trends <- vector()
      for(i in unique(years)) {
	    cat("  calculating moments for", i, "in", max(years), "\n")
		cat("    converting timeseries raster to matrix", "\n")
		  flush.console(); Sys.sleep(0.01)
		r <- as.data.frame(fcov[[which(years == i)]], cells = TRUE, na.rm = TRUE)[fcells,] 
          cells <- r[,1]
	        r <- r[,-1]
	    ctl.trends <- append(ctl.trends, median(as.matrix(r[which(cells %in% ctl.idx),]), na.rm = TRUE))  
		int.trends <- append(int.trends, median(as.matrix(r[which(cells %in% int.idx),]), na.rm = TRUE))
		out.trends <- append(out.trends, median(as.matrix(r[-which(cells %in% c(ctl.idx,int.idx)),]), na.rm = TRUE))
		remove(r, cells)
	  }
    plot(unique(years),int.trends, type="n", 
	     ylim=c(range(c(ctl.trends, int.trends, out.trends))),  
	     xlab = "year", ylab = "Median fCOV",  
         main = paste0(j, " fCOV  yearly trends"))
      lines(unique(years), ctl.trends, lty = 1, col = "red")
        points(unique(years), ctl.trends, pch=20, cex = 0.75, col = "red")
      lines(unique(years), int.trends, lty = 1, col = "green")
        points(unique(years), int.trends, pch=20, cex = 0.75, col = "green") 
      lines(unique(years), out.trends, lty = 3, col = "black")
        points(unique(years), out.trends, pch=20, cex = 0.75, col = "black")
	abline(lm(apply(rbind(ctl.trends, int.trends, out.trends), MARGIN=2, mean) ~ as.numeric(unique(years))), 
	       col="black", lty=6, lwd=0.5)		
    legend("bottomright", legend=c("Control", "Treatment", "Outside", "average trend"), lty=c(1,1,3, 6), 
           col=c("red", "green", "black", "black"))  
    remove(fcov, ctl.idx, int.idx, ctl.trends, int.trends, out.trends, fcov.con)	
  gc()	
  }		 
dev.off()


pdf(file.path("C:/evans/PFP/results", "LAI_yearly_trends2.pdf"), height=8.5, width=11)
  for(j in cty) {
    cat("Calculating trends for", j, "\n")
    setwd(file.path("C:/evans/PFP", j, "data"))
	  f <- rast("forest_250m.tif")
	    f[f < 1] <- NA
	      fcells <- as.data.frame(f, cells = TRUE, na.rm = TRUE)$cell 
      rn <- list.files(getwd(), "tif$") 
	    lai <- rast(grep("lai_250m_trend", rn, value = TRUE))
          dates <- as.Date(names(lai))
            years <- format(dates, "%Y")
    lai.con <- file.path("C:/evans/PFP", j, "model", "lai", "lai_250m_model_data.duckdb")
      con <- dbConnect(duckdb::duckdb(), lai.con, read_only=TRUE)
       db.table <- dbListTables(con)[1]
       dat <- dbGetQuery(con, paste0("SELECT * FROM ", db.table))
    dbDisconnect(con, shutdown = TRUE)
	  ctl.idx <- dat[which(dat$NAME == "control"),]$cell
      int.idx <- dat[which(dat$NAME != "control"),]$cell
	  remove(dat)
	ctl.trends <- vector()  
    int.trends <- vector()
    out.trends <- vector()
      for(i in unique(years)) {
	    cat("  calculating moments for", i, "in", max(years), "\n")
		cat("    converting timeseries raster to matrix", "\n")
		  flush.console(); Sys.sleep(0.01)
		r <- as.data.frame(lai[[which(years == i)]], cells = TRUE, na.rm = TRUE)[fcells,] 
          cells <- r[,1]
	        r <- r[,-1]
	    ctl.trends <- append(ctl.trends, median(as.matrix(r[which(cells %in% ctl.idx),]), na.rm = TRUE))  
		int.trends <- append(int.trends, median(as.matrix(r[which(cells %in% int.idx),]), na.rm = TRUE))
		out.trends <- append(out.trends, median(as.matrix(r[-which(cells %in% c(ctl.idx,int.idx)),]), na.rm = TRUE))
		remove(r, cells)
	  }
    plot(unique(years),int.trends, type="n", 
	     ylim=c(range(c(ctl.trends, int.trends, out.trends))),  
	     xlab = "year", ylab = "Median lai",  
         main = paste0(j, " lai  yearly trends"))
      lines(unique(years), ctl.trends, lty = 1, col = "red")
        points(unique(years), ctl.trends, pch=20, cex = 0.75, col = "red")
      lines(unique(years), int.trends, lty = 1, col = "green")
        points(unique(years), int.trends, pch=20, cex = 0.75, col = "green") 
      lines(unique(years), out.trends, lty = 3, col = "black")
        points(unique(years), out.trends, pch=20, cex = 0.75, col = "black")
	abline(lm(apply(rbind(ctl.trends, int.trends, out.trends), MARGIN=2, mean) ~ as.numeric(unique(years))), 
	       col="black", lty=6, lwd=0.5)		
    legend("bottomright", legend=c("Control", "Treatment", "Outside", "average trend"), lty=c(1,1,3, 6), 
           col=c("red", "green", "black", "black"))  
    remove(lai, ctl.idx, int.idx, ctl.trends, int.trends, out.trends, lai.con)	
  gc()	
  }		 
dev.off()

