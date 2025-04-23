suppressMessages(
  lapply(c("sf", "spatialEco", "terra",  "dplyr", 
           "duckdb", "DBI"), 
			require, character.only = TRUE))

cty <- c("brazil", "bhutan", "canada", "colombia", "costa_rica", "peru")

  switch(cty, 
    "brazil" = { mres = "250m" },  
    "bhutan" = { mres = "500m" },
	"canada" = { mres = "250m" }, 
	"colombia" = { mres = "250m" }, 
	"costa_rica" = { sub.months = c("12", "01", "02", "03", "04") }, 
	"peru" = { mres = "250m" } 
	) 
	
#*************************************************
#*************************************************
# Leaf Area Index (LAI)
#*************************************************
#*************************************************

for(j in cty) {
  cat("Calculating trends for", j, "\n")
  setwd(file.path("C:/evans/PFP", j, "model"))
    graph.dir = file.path("C:/evans/PFP/results", j, "timeseries_graphs") 
      dir.create(graph.dir, showWarnings = FALSE)

  pdf(file.path(graph.dir, "LAI_yearly_ecoregion_trends_new.pdf"), height=8.5, width=11)	  
  trends.out <- file.path(graph.dir, "lai_ecosystem_timeseries_new.csv")
  
  ts.db = file.path(getwd(), "lai", "database", "lai_250m_timeseries_data.duckdb") 
  dat.db = file.path(getwd(), "lai", "lai_250m_model_data.duckdb")
  
  con <- dbConnect(duckdb::duckdb(), dat.db, read_only=TRUE)
    db.table <- dbListTables(con)
      itype = as.data.frame(table(dbGetQuery(con, paste0("SELECT elu FROM ", db.table))[,1]))
    names(itype) <- c("ELU", "COUNT")
  dbDisconnect(con, shutdown = TRUE)
  
  itype <- itype[which(itype$COUNT >= 100),]  
  elist <- as.character(itype$ELU)

  rn <- list.files(file.path("C:/evans/PFP", j, "data"), "tif$", full.names=TRUE) 
    dates <- as.Date(names(rast(grep("lai_250m_trend", rn, value = TRUE))))
      years <- format(dates, "%Y")

  ctl.trends.list <- list()
  int.trends.list <- list()

  ct=0
  for(e in elist) {
    ct=ct+1
    db.query = paste0("SELECT * FROM data WHERE elu LIKE '%", e, "%'")
    con <- dbConnect(duckdb::duckdb(), dat.db, read_only=TRUE)
      db.table <- dbListTables(con)
      db.query <- gsub("data", db.table, db.query)
      dat <- dbGetQuery(con, db.query)
    dbDisconnect(con, shutdown = TRUE)
  
    ctl.idx <- dat[which(dat$NAME == "control"),]$cell
    int.idx <- dat[which(dat$NAME != "control"),]$cell
    remove(dat)

    ctl.df <- data.frame(matrix(, nrow = 1, ncol = 22))
	  names(ctl.df) <- 2000:2021
	int.df <- data.frame(matrix(, nrow = 1, ncol = 22))
	  names(int.df) <- 2000:2021	
	ctl.trends <- vector()
	int.trends <- vector()
	
    for(i in unique(years)) {
	  cat(j, " - calculating median for", i, "in", max(years), "(", ct, "of", length(elist), "ecosystems )",  "\n")
	    flush.console(); Sys.sleep(0.01)	
	  y <- dates[which(years == i)]
	    y <- y[grep(paste(sub.months, collapse="|"), y)]  
	      y <- paste0("D", gsub("-", "", as.character(y)))
	  con <- dbConnect(duckdb::duckdb(dbdir = ts.db), read_only = FALSE)
	    db.table <- dbListTables(con)
		tdat <- dbGetQuery(con, paste0("SELECT COLUMNS('(", paste(c("cell", y), collapse ="|"), ")') FROM ", db.table))
        # dbGetQuery(con, paste0("SELECT * FROM ", db.table, " WHERE rowid BETWEEN ", 1, " AND " 2))   		
      dbDisconnect(con, shutdown = TRUE)
	  cells <- tdat$cell
	    tdat <- tdat[,-1]
	  ctl.trends <- append(ctl.trends, median(as.matrix(tdat[which(cells %in% ctl.idx),]), na.rm = TRUE))  
	  int.trends <- append(int.trends, median(as.matrix(tdat[which(cells %in% int.idx),]), na.rm = TRUE))
	    remove(tdat, cells, y)
      gc()	  
	}
	ctl.trends.list[[e]] <- ctl.trends 
	int.trends.list[[e]] <- int.trends 
    ctl.df[1,] <- ctl.trends 
	  ctl.df <- data.frame(country=j, ecosystem=e, metric = "fcov", intervention="control", ctl.df)
	int.out <- file.path(graph.dir, "fcov_treatment_ecosystem_timeseries.csv")
    int.df[1,] <- int.trends 
	  int.df <- data.frame(country=j, ecosystem=e, metric = "fcov", intervention="treatment", int.df)
	trends.df <- rbind(int.df, ctl.df)  
    write.table(trends.df, trends.out, sep = ",", row.names = FALSE, append = TRUE, 
                col.names = !file.exists(trends.out))
    remove(ctl.df, int.df, ctl.idx, int.idx, trends.df)				
  gc()
  }
    plot(unique(years), ctl.trends.list[[1]], type="n", 
	     ylim = range(c(unlist(ctl.trends.list), unlist(int.trends.list))),  
	     xlab = "year", ylab = "Median LAI by Ecosystem Type",  
         main = paste0(j, " LAI  yearly trends"))
      
	  for(i in 1:length(ctl.trends.list)){
	    lines(unique(years), ctl.trends.list[[i]], lty=3, col="red")
          points(unique(years), ctl.trends.list[[i]], pch=20, cex = 0.75, col="red")
	  }
	  for(i in 1:length(int.trends.list)){
	    lines(unique(years), int.trends.list[[i]], col="darkgreen")
          points(unique(years), int.trends.list[[i]], pch=20, cex = 0.75, col="darkgreen")
	  }	  
	legend("bottomright", legend=c("treatment", "control"), lty=c(1,3), 
	       pch=c(20,20), col=c("darkgreen", "red"))  
  dev.off()
}

lapply(int.trends.list, kendall)


#*************************************************
#*************************************************
# Fractional Cover (fCOV)
#*************************************************
#*************************************************

for(j in cty) {
  cat("Calculating trends for", j, "\n")
  setwd(file.path("C:/evans/PFP", j, "model"))
    graph.dir = file.path("C:/evans/PFP/results", j, "timeseries_graphs") 
      dir.create(graph.dir, showWarnings = FALSE)

  pdf(file.path(graph.dir, "fCOV_yearly_ecoregion_trends.pdf"), height=8.5, width=11)	  
	  
  ts.db = file.path(getwd(), "fcov", "database", "fcov_500m_timeseries_data.duckdb") 
  dat.db = file.path(getwd(), "fcov", "fcov_500m_model_data.duckdb")
  trends.out <- file.path(graph.dir, "fcov_ecosystem_timeseries.csv")
  con <- dbConnect(duckdb::duckdb(), dat.db, read_only=TRUE)
    db.table <- dbListTables(con)
      itype = as.data.frame(table(dbGetQuery(con, paste0("SELECT elu FROM ", db.table))[,1]))
    names(itype) <- c("ELU", "COUNT")
  dbDisconnect(con, shutdown = TRUE)
  
  itype <- itype[which(itype$COUNT >= 100),]  
  elist <- as.character(itype$ELU)

  rn <- list.files(file.path("C:/evans/PFP", j, "data"), "tif$", full.names=TRUE) 
    dates <- as.Date(names(rast(grep("fcov_500m_trend", rn, value = TRUE))))
      years <- format(dates, "%Y")

  ctl.trends.list <- list()
  int.trends.list <- list()

  ct=0
  for(e in elist) {
    ct=ct+1
    db.query = paste0("SELECT * FROM data WHERE elu LIKE '%", e, "%'")
    con <- dbConnect(duckdb::duckdb(), dat.db, read_only=TRUE)
      db.table <- dbListTables(con)
      db.query <- gsub("data", db.table, db.query)
      dat <- dbGetQuery(con, db.query)
    dbDisconnect(con, shutdown = TRUE)
  
    ctl.idx <- dat[which(dat$NAME == "control"),]$cell
    int.idx <- dat[which(dat$NAME != "control"),]$cell
    remove(dat)

    ctl.df <- data.frame(matrix(, nrow = 1, ncol = 22))
	  names(ctl.df) <- 2000:2021
	int.df <- data.frame(matrix(, nrow = 1, ncol = 22))
	  names(int.df) <- 2000:2021	
	ctl.trends <- vector()
	int.trends <- vector()
	
    for(i in unique(years)) {
	  cat(j, " - calculating median for", i, "in", max(years), "(", ct, "of", length(elist), "ecosystems )",  "\n")
	    flush.console(); Sys.sleep(0.01)	
	  y <- dates[which(years == i)]
	    y <- y[grep(paste(sub.months, collapse="|"), y)]
	      y <- paste0("D", gsub("-", "", as.character(y)))
	  con <- dbConnect(duckdb::duckdb(dbdir = ts.db), read_only = FALSE)
	    db.table <- dbListTables(con)
		tdat <- dbGetQuery(con, paste0("SELECT COLUMNS('(", paste(c("cell", y), collapse ="|"), ")') FROM ", db.table))  
      dbDisconnect(con, shutdown = TRUE)
	  cells <- tdat$cell
	    tdat <- tdat[,-1]
	  ctl.trends <- append(ctl.trends, median(as.matrix(tdat[which(cells %in% ctl.idx),]), na.rm = TRUE))  
	  int.trends <- append(int.trends, median(as.matrix(tdat[which(cells %in% int.idx),]), na.rm = TRUE))
	    remove(tdat, cells, y)
      gc()	  
	}
	ctl.trends.list[[e]] <- ctl.trends 
	int.trends.list[[e]] <- int.trends 
    ctl.df[1,] <- ctl.trends 
	  ctl.df <- data.frame(country=j, ecosystem=e, metric = "fcov", intervention="control", ctl.df)
	int.out <- file.path(graph.dir, "fcov_treatment_ecosystem_timeseries.csv")
    int.df[1,] <- int.trends 
	  int.df <- data.frame(country=j, ecosystem=e, metric = "fcov", intervention="treatment", int.df)
	trends.df <- rbind(int.df, ctl.df)  
    write.table(trends.df, trends.out, sep = ",", row.names = FALSE, append = TRUE, 
                col.names = !file.exists(trends.out))
    remove(ctl.df, int.df, ctl.idx, int.idx, trends.df)				
  gc()
  }
    plot(unique(years), ctl.trends.list[[1]], type="n", 
	     ylim = range(c(unlist(ctl.trends.list), unlist(int.trends.list))),  
	     xlab = "year", ylab = "Median fCOV by Ecosystem Type",  
         main = paste0(j, " fCOV  yearly trends"))
      
	  for(i in 1:length(ctl.trends.list)){
	    lines(unique(years), ctl.trends.list[[i]], lty=3, col="red")
          points(unique(years), ctl.trends.list[[i]], pch=20, cex = 0.75, col="red")
	  }
	  for(i in 1:length(int.trends.list)){
	    lines(unique(years), int.trends.list[[i]], col="darkgreen")
          points(unique(years), int.trends.list[[i]], pch=20, cex = 0.75, col="darkgreen")
	  }	  
	legend("bottomright", legend=c("treatment", "control"), lty=c(1,3), 
	       pch=c(20,20), col=c("darkgreen", "red"))  
  dev.off()
}

lapply(int.trends.list, kendall)



#*************************************************
#*************************************************
# Aggregrate rasters by year
#*************************************************
#*************************************************
suppressMessages(
  lapply(c("sf", "spatialEco", "terra",  "dplyr", "rts",  
           "duckdb", "DBI"), 
			require, character.only = TRUE))

cty <- c("brazil", "bhutan", "canada", "colombia", "costa_rica", "peru")[-1]
  for(j in cty) {
    root <- file.path("C:/evans/PFP", j)
    setwd(file.path(root, "data"))
      out.dir <- file.path(root, "data", "yearly_data")
        dir.create(out.dir, showWarnings = FALSE)
    rn <- list.files(getwd(), "tif$", full.names=TRUE) 
    
    if(!file.exists(file.path(out.dir, "fcov_500m_yearly.tif"))) {
      cat("Calculating yearly fcov rasters for", j, "\n")
        flush.console(); Sys.sleep(0.01)	
      fcov <- rast(grep("fcov_500m_trend", rn, value = TRUE))
        if(j == "peru") {
          dates <- as.Date(names(rast(grep("lai_250m_trend", rn, value = TRUE))))[-1]
    	    names(fcov) <- dates
        } else {
          dates <- as.Date(names(fcov))	
    	}
      fcov <- rts(fcov, dates)
      system.time({ fcov.yr <- apply.yearly(fcov, median)@raster })
        names(fcov.yr) <- 2000:2021
      writeRaster(fcov.yr, file.path(out.dir, "fcov_500m_yearly.tif"))
    }
    if(!file.exists(file.path(out.dir, "lai_250m_yearly.tif"))) {
      cat("Calculating yearly lai rasters for", j, "\n")
        flush.console(); Sys.sleep(0.01)	  
      lai <- rast(grep("lai_250m_trend", rn, value = TRUE))
	    dates <- as.Date(names(lai))
          lai <- rts(lai, dates) 
      system.time({ lai.yr <- apply.yearly(lai, median)@raster })
        names(lai.yr) <- 2000:2021
      writeRaster(lai.yr, file.path(out.dir, "lai_250m_yearly.tif")) 
    }
    tmpFiles(current=TRUE, orphan=TRUE, old=TRUE, remove=TRUE)
  gc()	
  }
