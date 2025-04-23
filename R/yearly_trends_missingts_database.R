suppressMessages(
  lapply(c("sf", "spatialEco", "terra",  "dplyr", 
           "duckdb", "DBI"), 
			require, character.only = TRUE))

cty <- c("brazil", "bhutan", "canada", "colombia", "costa_rica", "peru")[6]

#*************************************************
#*************************************************
# Leaf Area Index (LAI)
#*************************************************
#*************************************************

j = cty
  cat("Calculating trends for", j, "\n")
  setwd(file.path("C:/evans/PFP", j, "model"))
    graph.dir = file.path("C:/evans/PFP/results", j, "timeseries_graphs") 
      dir.create(graph.dir, showWarnings = FALSE)
  pdf(file.path(graph.dir, "lai_yearly_ecoregion_trends_peru.pdf"), height=8.5, width=11)	  
  dat.db = file.path(getwd(), "lai", "lai_250m_model_data.duckdb")
  trends.out <- file.path(graph.dir, "lai_ecosystem_timeseries.csv")
  con <- dbConnect(duckdb::duckdb(), dat.db, read_only=TRUE)
    db.table <- dbListTables(con)
      itype = as.data.frame(table(dbGetQuery(con, paste0("SELECT elu FROM ", db.table))[,1]))
    names(itype) <- c("ELU", "COUNT")
  dbDisconnect(con, shutdown = TRUE)
    itype <- itype[which(itype$COUNT >= 100),]  
    elist <- as.character(itype$ELU)
  rn <- list.files(file.path("C:/evans/PFP", j, "data"), "tif$", full.names=TRUE)
    lai <- rast(grep("lai_250m_trend", rn, value = TRUE))  
      dates <- as.Date(names(lai))	
        years <- format(dates, "%Y")
  ct=0
  for(e in elist) {
    ct=ct+1
    db.query = paste0("SELECT * FROM data WHERE elu LIKE '%", e, "%'")
    con <- dbConnect(duckdb::duckdb(), dat.db, read_only=TRUE)
      db.table <- dbListTables(con)
      db.query <- gsub("data", db.table, db.query)
      dat <- dbGetQuery(con, db.query)
    dbDisconnect(con, shutdown = TRUE)
	ctl <- dat[which(dat$NAME == "control"),][c("X", "Y")]
	treat <- dat[which(dat$NAME != "control"),][c("X", "Y")]
    remove(dat)
    ctl.trends.list <- list()
    int.trends.list <- list()
    ctl.df <- data.frame(matrix(, nrow = 1, ncol = 22))
	  names(ctl.df) <- 2000:2021
	int.df <- data.frame(matrix(, nrow = 1, ncol = 22))
	  names(int.df) <- 2000:2021	
	ctl.trends <- vector()
	int.trends <- vector()
      for(i in unique(years)) {
	    cat(j, " - calculating median for", i, "in", max(years), "(", ct, "of", length(elist), "ecosystems )",  "\n")
	      flush.console(); Sys.sleep(0.01)	
	    ctl.trends <- append(ctl.trends, median(as.matrix(extract(lai[[which(years == i)]], ctl[,c("X","Y")])[,-1]), na.rm = TRUE))  
	    int.trends <- append(int.trends, median(as.matrix(extract(lai[[which(years == i)]], treat[,c("X","Y")])[,-1]), na.rm = TRUE))
        gc()	  
	  }
	ctl.trends.list[[e]] <- ctl.trends 
	int.trends.list[[e]] <- int.trends 
    ctl.df[1,] <- ctl.trends 
	  ctl.df <- data.frame(country=j, ecosystem=e, metric = "lai", intervention="control", ctl.df)
	int.out <- file.path(graph.dir, "lai_treatment_ecosystem_timeseries.csv")
    int.df[1,] <- int.trends 
	  int.df <- data.frame(country=j, ecosystem=e, metric = "lai", intervention="treatment", int.df)
	trends.df <- rbind(int.df, ctl.df)  
    write.table(trends.df, trends.out, sep = ",", row.names = FALSE, append = TRUE, 
                col.names = !file.exists(trends.out))
    remove(ctl.df, int.df, ctl, treat, trends.df)				
  gc()
  }
    plot(unique(years), ctl.trends.list[[1]], type="n", 
	     ylim = range(c(unlist(ctl.trends.list), unlist(int.trends.list))),  
	     xlab = "year", ylab = "Median lai by Ecosystem Type",  
         main = paste0(j, " lai  yearly trends"))
      
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


tdat <- read.csv(trends.out)
  treat <- tdat[tdat$intervention == "treatment",][-c(1:4)]
  ctl <- tdat[tdat$intervention == "control",][-c(1:4)]
  pdf(file.path(graph.dir, "lai_yearly_ecoregion_trends_peru.pdf"), height=8.5, width=11)	  
    plot(unique(years), treat[1,], type="n", 
	     ylim = range(rbind(treat, ctl)),  
	     xlab = "year", ylab = "Median lai by Ecosystem Type",  
         main = paste0(j, " lai  yearly trends"))
      
	  for(i in 1:nrow(ctl)){
	    lines(unique(years), ctl[i,], lty=3, col="red")
          points(unique(years), ctl[i,], pch=20, cex = 0.75, col="red")
	  }
	  for(i in 1:nrow(treat)){
	    lines(unique(years), treat[i,], col="darkgreen")
          points(unique(years), treat[i,], pch=20, cex = 0.75, col="darkgreen")
	  }	  
	legend("bottomright", legend=c("treatment", "control"), lty=c(1,3), 
	       pch=c(20,20), col=c("darkgreen", "red"))  
dev.off()


#*************************************************
#*************************************************
# Fractional Cover (fCOV)
#*************************************************
#*************************************************

j = cty
  cat("Calculating trends for", j, "\n")
  setwd(file.path("C:/evans/PFP", j, "model"))
    graph.dir = file.path("C:/evans/PFP/results", j, "timeseries_graphs") 
      dir.create(graph.dir, showWarnings = FALSE)
  pdf(file.path(graph.dir, "fCOV_yearly_ecoregion_trends.pdf"), height=8.5, width=11)	  
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
    fcov <- rast(grep("fcov_500m_trend", rn, value = TRUE))  
      if(j == "peru") {
        dates <- as.Date(names(rast(grep("lai_250m_trend", rn, value = TRUE))))[-1]
          names(fcov) <- dates
      } else {
        dates <- as.Date(names(fcov))	
      }
      years <- format(dates, "%Y")
  ct=0
  for(e in elist) {
    ct=ct+1
    db.query = paste0("SELECT * FROM data WHERE elu LIKE '%", e, "%'")
    con <- dbConnect(duckdb::duckdb(), dat.db, read_only=TRUE)
      db.table <- dbListTables(con)
      db.query <- gsub("data", db.table, db.query)
      dat <- dbGetQuery(con, db.query)
    dbDisconnect(con, shutdown = TRUE)
	ctl <- dat[which(dat$NAME == "control"),][c("X", "Y")]
	treat <- dat[which(dat$NAME != "control"),][c("X", "Y")]
    remove(dat)
    ctl.trends.list <- list()
    int.trends.list <- list()
    ctl.df <- data.frame(matrix(, nrow = 1, ncol = 22))
	  names(ctl.df) <- 2000:2021
	int.df <- data.frame(matrix(, nrow = 1, ncol = 22))
	  names(int.df) <- 2000:2021	
	ctl.trends <- vector()
	int.trends <- vector()
      for(i in unique(years)) {
	    cat(j, " - calculating median for", i, "in", max(years), "(", ct, "of", length(elist), "ecosystems )",  "\n")
	      flush.console(); Sys.sleep(0.01)	
	    ctl.trends <- append(ctl.trends, median(as.matrix(extract(fcov[[which(years == i)]], ctl[,c("X","Y")])[,-1]), na.rm = TRUE))  
	    int.trends <- append(int.trends, median(as.matrix(extract(fcov[[which(years == i)]], treat[,c("X","Y")])[,-1]), na.rm = TRUE))
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
    remove(ctl.df, int.df, ctl, treat, trends.df)				
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


trends.out <- file.path(graph.dir, "fcov_ecosystem_timeseries.csv")
tdat <- read.csv(trends.out)
  treat <- tdat[tdat$intervention == "treatment",][-c(1:4)]
  ctl <- tdat[tdat$intervention == "control",][-c(1:4)]
  pdf(file.path(graph.dir, "fcov_yearly_ecoregion_trends_peru.pdf"), height=8.5, width=11)	  
    plot(unique(years), treat[1,], type="n", 
	     ylim = range(rbind(treat, ctl)),  
	     xlab = "year", ylab = "Median fcov by Ecosystem Type",  
         main = paste0(j, " fcov  yearly trends"))
      
	  for(i in 1:nrow(ctl)){
	    lines(unique(years), ctl[i,], lty=3, col="red")
          points(unique(years), ctl[i,], pch=20, cex = 0.75, col="red")
	  }
	  for(i in 1:nrow(treat)){
	    lines(unique(years), treat[i,], col="darkgreen")
          points(unique(years), treat[i,], pch=20, cex = 0.75, col="darkgreen")
	  }	  
	legend("bottomright", legend=c("treatment", "control"), lty=c(1,3), 
	       pch=c(20,20), col=c("darkgreen", "red"))  
dev.off()