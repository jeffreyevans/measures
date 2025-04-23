suppressMessages( lapply(c("terra",  "duckdb", "DBI"), 
                  require, character.only = TRUE))
			
#***************************
# Set working directory
path = "C:/evans/PFP"
country <- c("brazil", "bhutan", "canada", "colombia", "costa_rica", "peru")[6]

root = file.path(path, country)
setwd(file.path(root, "data"))
  mdl.dir = file.path(root, "model")
  dat.dir = file.path(root, "data")

#***************************
# Metric and resoultion; 
metric.type <- c("lai", "fcov")[2]  
  switch(metric.type, 
    "lai" = { mres = "250m" },  
    "fcov" = { mres = "500m" }) 
  switch(metric.type, 
    "lai" = { diff.thres = 0.5 },  
    "fcov" = { diff.thres = 0.01 }) 

#*************************************************
# set database connections
con.dat <- file.path(mdl.dir, metric.type,paste0(metric.type, "_",mres, "_model_data.duckdb"))        # final model data

#*************************************************
# read metric raster

if(metric.type == "lai"){
  metric <- rast(file.path(dat.dir, paste0(toupper(metric.type), "_", mres, "_trend.tif")))
    metric.max <- max(global(metric, max, na.rm=TRUE)[,1])
} else {
  metric.max <- 1
}

#*************************************************
# read model data
con <- dbConnect(duckdb::duckdb(), con.dat, read_only=TRUE)
  db.table <- dbListTables(con)
  dat <- dbGetQuery(con, paste0("SELECT * FROM ", db.table))
dbDisconnect(con, shutdown = TRUE)

if(country == "canada") {
dat$NAME <- gsub("'","", dat$NAME)
  dat$NAME <- gsub("/", " ", dat$NAME, fixed = TRUE)
    dat$NAME <- gsub("?", " ", dat$NAME, fixed = TRUE)
	  dat$NAME <- gsub("&", "and", dat$NAME, fixed = TRUE)
}

#*************************************************
# Calculate Gain and Loss
current <- dat$current 
baseline <- dat$pre 

gain <- current - baseline
  gain[gain == 0] <- 0.01
    gain[gain < 0] <- 0
loss <- current - baseline
  loss[loss == 0] <- -0.01
    loss[loss > 0] <- 0
      loss <- loss * -1    

cpct <- (current / baseline) - 1
  lt.idx <- which(cpct < 0)
    for(i in lt.idx) {
      if( abs(current[i] - baseline[i]) <= diff.thres ) {
        cpct[i] <- 0
      } else {
        cpct[i] <- (current[i] / baseline[i]) / metric.max	  
	  }
    }	
  gt.idx <- which(cpct > 1)
    for(i in gt.idx) {
      if( abs(current[i] - baseline[i]) <= diff.thres ) {
        cpct[i] <- 0
      } else {
        cpct[i] <- (current[i] - baseline[i]) / metric.max	  
	  }
    }

summary(gain)
summary(loss)
summary(cpct)

#gt.idx <- which(cpct > 1)
#  if(length(gt.idx) > 1) 
#    current[gt.idx] / baseline[gt.idx) * 0.1	
#summary(cpct)	
#  plot(density(cpct))
#

dat$gain <- gain
dat$loss <- loss
dat$cpct <- cpct

#*************************************************
# Write new results
con <- dbConnect(duckdb::duckdb(), con.dat, read_only=FALSE)
  duckdb::dbWriteTable(con, paste0(db.table), dat, overwrite = TRUE)
dbDisconnect(con, shutdown = TRUE)   	

remove(dat, cpct, gain, loss, current, baseline)

#head(dat[gt.idx,][,c("pre", "current")])
