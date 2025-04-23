#*************************************************
# Causal Ensemble and Double Debiased ML Intervention Impact Model
#
# Expected input data
#   * "country.gpkg - interventions" are vector polygons of the
#      protected areas and PFP's (if applicable). These polygons are 
#      also used to specify the CATE (mean treatment effect) estimates.    
#  * "model_data.duckdb" - Duck Database file 
#     Observations of treatments and controls (pixel centers). The "cell" column corresponds
#      to the cell index in the mask raster. The "DESIG" column indicates intervention type as
#      well as "control" observations. The "NAME" column is the name of the intervention and used 
#      to specified models. Data also includes the dependent and independet variables. The "source_cell" 
#      is what relates controls back to interventions.   
#        
#      This requires output from the BuildModelData.R script 
# 
# Model outputs
#   * "y_mres_effect_size.tif" and the hetrogenious effect size (pixel-level) estimates
#      written as a 32-bit (FLT4S) floating point raster with "effect.size", "lci", "uci"
#      layers representing effect sizes and upper/lower confidence bounds
#
#*************************************************
# Required libraries
suppressMessages(
  lapply(c("sf", "spatialEco", "terra",  "dplyr", "fst", "duckdb", 
           "DBI", "grf", "rts", "ggplot2", "caret", "ranger", "ggplot2",   
		   "rfUtilities", "stars", "data.table", "doParallel"), 
			require, character.only = TRUE))
			
path = "C:/evans/PFP"
country <- c("brazil", "bhutan", "canada", "colombia", "costa_rica", "peru")[5]

#*************************************************
#*************************************************
# Model arguments

#***************************
# Metric and resoultion; 
#   lai 300m 2014-2024, lai 250m 2000 - 2021
#   fcov 300m 2014-2024, fcov 500m 2000 - 2021
metric.type <- c("lai", "fcov")[1]  # Metric
  switch(metric.type, 
    "lai" = { mres = "250m" },  
    "fcov" = { mres = "500m" }) 

#***************************
# Dependent variable;
#   current - current condition
#   cpct - percent that is current
#   tsa - temporal volume (area under the curve)
y = c("current", "tsa", "gain", "loss", "cpct")[4]         # Dependent variable
  yname = paste0(metric.type, "_", y)

#***************************
# paramters
group.by = c("pa", "elu")[1]               # What to aggregrate models by, protected areas or ELU's
nboot = 1001                               # Number of HTE Bootstrap replicates
p.sig = c(TRUE, FALSE)[1]                  # run parameter signifiance test
rm.cor = c(TRUE, FALSE)[1]                 # Evaluate and remove multi-colinear parameters
nt = round(parallel::detectCores()-1)      # sets number of cores for multiprocessing
mns = 10                                   # Minimum node size
max.size = 200000                          # Size of a model before triggering subsampling
mdl.size = 20000                           # size of subsampled models
write.db = c(TRUE, FALSE)[2]               # Write estimates to DuckDB database
table.idx = c(1,2)[1]                      # For Brazil LAI, which table 1-mountains, 2-plains
  if(country != "brazil") table.idx = 1 

# Names of dependet variables
dep.vars = c("current", "cpct", "tsa", "gain", "loss", "diff", "pchg")
		  
#*************************************************
#*************************************************
# Set working environment 
root = file.path(path, country)
setwd(file.path(root, "data"))
  mdl.dir = file.path(root, "model")
  dat.dir = file.path(root, "data")
  pred.dir = file.path(mdl.dir, metric.type, y)
    dir.create(pred.dir, showWarnings = FALSE)

#*************************************************
#*************************************************
# Functions
group.data.frame <- function(x, n = 10000, randomize = TRUE, mss = NULL, 
                             reallocate = c("fill", "back"), min.bin.size = NULL) {
  if(is.null(min.bin.size)) bin.bin.size = n
  if(inherits(x, "data.frame")) {  
    xf <- 1:nrow(x)
  } else {
    stop("x must be a data.frame object")
  }  
  if(randomize){ 
    xf <- sample(xf, length(xf))
  }
  s <- split(x, ceiling(seq_along(xf) / n))
    cts <- unlist(lapply(s, nrow))
  if(reallocate == "fill" & min(cts) < n) {
    if(min(cts) >= min.bin.size) {
    target.n <- diff(c(min(cts),n))
    s[[which.min(cts)]] <- rbind(s[[which.min(cts)]], 
	  do.call(rbind, lapply(1:target.n, \(i) {
        s[[sample(1:(length(s)-1), 1)]][sample(1:n, 1),] 
      })))
	} else {
      s <- s[-which.min(cts)]
    }	
  }
  if(!is.null(mss) & min(cts) < mss & min(cts) < n){	
    if(reallocate == "back" & min(cts) < n) {
      mdf <- s[[which.min(cts)]]
	    s <- s[-which.min(cts)] 
        for(r in 1:nrow(mdf)){
          rdf <- sample(length(s), 1)
          s[[rdf]] <- rbind(s[[rdf]], mdf[r,]) 
	    }
	  }
	}
  return(s) 
}

#*************************************************
# set database connections
con.dat <- file.path(mdl.dir, metric.type,paste0(metric.type, "_",mres, "_model_data.duckdb"))        # final model data

#*************************************************
#*************************************************
# read interventions, pa and pfp vectors 
#  ( lyrs <- st_layers(paste0(country, ".gpkg"))$name )
bdy <- st_read(paste0(country, ".gpkg"), "boundary" )
  if(country %in% c("brazil", "colombia", "peru")) {
    pfp <- st_cast(st_read(paste0(country, ".gpkg"), "PFP"), "POLYGON")
        st_geometry(pfp) <- "geometry" 	
    pa <- st_cast(st_read(paste0(country, ".gpkg"), "interventions"), "POLYGON")
        st_geometry(pa) <- "geometry"  
  } else {
    pa <- st_cast(st_read(paste0(country, ".gpkg"), "protected_areas"), "POLYGON")
      st_geometry(pa) <- "geometry"
  }
#
#*************************************************
#*************************************************
# read and que model data (locations/interventions, dependent, independent)  
##  database connections
#     con.dat - Final combined data "model_data.duckdb"
# Read attributed model data

# Selection intervention units with >= 100 observsations
if(group.by == "pa") {  
  con <- dbConnect(duckdb::duckdb(), con.dat, read_only=TRUE)
    db.table <- dbListTables(con)[table.idx]  
    n = dbGetQuery(con, paste0("SELECT COUNT(*) FROM ", db.table))[,1,1] 
    itype = as.data.frame(table(dbGetQuery(con, paste0("SELECT NAME FROM ", db.table))[,1]))
	  names(itype) <- c("NAME", "COUNT")
  dbDisconnect(con, shutdown = TRUE)
    itype <- itype[which(itype$COUNT >= 100),]  
      itype <- itype[order(itype$COUNT),]
        ulist <- as.character(itype$NAME)
	      ulist <- ulist[-which(itype$NAME == "control")]
	cat("Model sizes;", "\n")
	  print(itype)  			  
} else if(group.by == "elu"){
  poly_es = FALSE
  con <- dbConnect(duckdb::duckdb(), con.dat, read_only=TRUE)
    db.table <- dbListTables(con)[table.idx]  
    itype = as.data.frame(table(dbGetQuery(con, paste0("SELECT elu FROM ", db.table))[,1]))
	  names(itype) <- c("ELU", "COUNT")
  dbDisconnect(con, shutdown = TRUE)
    itype <- itype[which(itype$COUNT >= 100),]  
      itype <- itype[order(itype$COUNT),]
        ulist <- as.character(itype$ELU)
	cat("Model sizes;", "\n")
	  print(itype)
}

#*************************************************
# Create polygon ATE table and read reference raster 
# for HTE estimates, preds correspond to cell ids
cf.ate <- list()

es.name = file.path(pred.dir, paste(y, mres, "effect_size.tif", sep="_"))
if(!file.exists(es.name)) {
  es.raster <- rast(paste0("forest_",mres,".tif"))
    es.raster <- rep(es.raster, 5)
      es.raster[] <- NA
        names(es.raster) <- c("effect.size", "variance", "unscaled.es", "LCI", "UCI")    
  cat("**** Writing HTE raster to", es.name, "\n")
  es.raster <- writeRaster(es.raster, es.name, datatype="FLT4S", overwrite=TRUE)
} else {
  es.raster <- rast(es.name)
    #es.raster <- es.raster[[which(names(es.raster) %in% c("effect.size", "variance", "unscaled.es", "LCI", "UCI"))]]
}

#*************************************************
# Read controls
db.query = paste0("SELECT * FROM data WHERE NAME LIKE '%", "control", "%'")
con <- dbConnect(duckdb::duckdb(), con.dat, read_only=TRUE)
  db.table <- dbListTables(con)[table.idx]
  db.query <- gsub("data", db.table, db.query)
  controls <- dbGetQuery(con, db.query)
dbDisconnect(con, shutdown = TRUE)

drop.idx <- which(names(controls) %in% c("NAME_IDX", "patch")) 
  if(length(drop.idx) > 0) controls <- controls[-drop.idx,]

controls <- controls[,-which(names(controls) %in% c("X", "Y"))]
  names(controls)[which(names(controls) %in% y)] <- "y"
    controls <- controls[,-which(names(controls) %in% dep.vars)]

parms <- names(controls)[which(names(controls) %in% "pre"):ncol(controls)]
  parms <- parms[-which(parms == "y")]

#******************************************************************
#******************************************************************
# **** Start of for loop for processing each intervention type **** 
# ulist <- ulist[ct:length(ulist)]
# ulist = ulist[c(210,212,211)]  # colombia missing
# ulist = ulist = ulist[52]      # peru missing

ct = 0
  for(j in ulist) {
    ct = ct+1
    cat("***************************************************", "\n")
    cat("**** Causal Impact Ensemble Model for:", j, "in", country, "\n")
	cat("****   Country", country, "\n")
	cat("****   Intervention:", j, "\n")
    cat("****   Aggregrating models by:", group.by, "\n")	
	cat("****   Model iteration", ct, "in", length(ulist), "\n")	
    cat("****   Metric:", metric.type, "\n")
    cat("****   Dependent variable is:", y, "\n")
    cat("\n")	

    #################################
    #################################  
    # Build model data from SQL query
	switch(group.by, 
	  "pa" = { db.query = paste0("SELECT * FROM data WHERE NAME LIKE '%", j, "%'")},
	  "elu" = {db.query = paste0("SELECT * FROM data WHERE elu LIKE '%", j, "%'")})
	con <- dbConnect(duckdb::duckdb(), con.dat, read_only=TRUE)
	  db.table <- dbListTables(con)[table.idx]
	  db.query <- gsub("data", db.table, db.query)
	  isub <- dbGetQuery(con, db.query)
	dbDisconnect(con, shutdown = TRUE)
  	
	drop.idx <- which(names(isub) %in% c("NAME_IDX", "patch")) 
      if(length(drop.idx) > 0) isub <- isub[-drop.idx,]
	ctl.idx <- which(isub$NAME %in% "control")
      if(length(ctl.idx) > 0) isub <- isub[-ctl.idx,]
	if(!exists("isub") | nrow(isub) < 100) { next }
	isub$treatment <- 1
      names(isub)[which(names(isub) %in% y)] <- "y"
        isub <- isub[,-which(names(isub) %in% dep.vars)]
	      didx <- which(duplicated(isub$cell)) 
            if(length(didx) > 0) isub <- isub[-didx,] 
	if(nrow(isub)*2 > max.size) {
      cat("****  Randomly partitioning into multiple models", "\n")
	    flush.console(); Sys.sleep(0.01)		  
	    mdata <- group.data.frame(isub, n = mdl.size/2, mss = 100, reallocate = "fill", min.bin.size = mdl.size/2*0.25)
    	  mdata <- lapply(mdata, \(m) {
			#names(m)[which(names(m) %in% y)] <- "y"
           	  cs <- controls[which(m$cell %in% controls$source_cell),]
                cs$treatment <- 0
              return(rbind(m[,c("cell", "y", "treatment", parms)],
    		         cs[,c("cell", "y", "treatment", parms)]))
    		})
	  remove(isub)
	  cat("****    partitioning resulted in", length(mdata), "models", "\n")
	  cat("\n")  
    } else {
    csub <- controls[which(isub$cell %in% controls$source_cell),]
	 csub$treatment <- 0
	    if(nrow(csub) < 20) {
	      warning(j, "does not have enough observations, skipping model")
	        next	    
		}
	    if(nrow(csub) > (nrow(isub) + (nrow(isub) * 0.10))) {  
	      csub <- csub[sample(1:nrow(csub), nrow(isub)),]
	    }
        mdata <- rbind(isub[,c("cell", "y", "treatment", parms)], 
		               csub[,c("cell", "y", "treatment", parms)])	
      remove(csub,isub)
	}
	  mdata <- na.omit(mdata)
	gc(verbose = FALSE, full = TRUE)
	
    #################################
    #################################  
    # Screening for collinear and multi-collinear variables
	if(rm.cor | p.sig) {
	  if(inherits(mdata, "list")){
	    X <- do.call(rbind, mdata)
	  } else {
	    X <- mdata
	  }
	  rm.idx <- which(names(X) %in% c("cell", "y")) 
	    if(length(rm.idx) > 0) X <- X[,-rm.idx]
	      X <- na.omit(X)	  
	}
	
    if(rm.cor) {
	options(warn=-1)
	cat("****  Screening for collinear and multi-collinear variables", "\n")
	  flush.console(); Sys.sleep(0.01)		
      all.vars <- names(X)[-1]
        cl.vars <- try(suppressMessages(collinear(X[,all.vars], p = 0.85,  
  	                   nonlinear = FALSE, p.value = 0.001) ) )				 
  	  if(exists("cl.vars")){
        if(length(cl.vars) > 0)
  		  if("trend" %in% cl.vars) cl.vars <- cl.vars[-which(cl.vars %in% "trend")] 
            X <- X[,-which(names(X) %in% cl.vars)]
  	  }	
      mc <- suppressMessages(multi.collinear(X[,-1], p = 0.05, perm = TRUE, n = 99))
        mc.vars <- mc[which( mc$frequency > 5 ),]$variables 
          if(length(mc.vars) > 0){
            if("trend" %in% mc.vars) mc.vars <- mc.vars[-which(mc.vars %in% "trend")] 
              X <- X[,-which(names(X) %in% mc.vars)]
  		}
	  cat("\n")
      cat("The following variables were dropped:",
        all.vars[which(is.na(match(all.vars, names(X))))], "\n")
      cat("The retained variables are:", 
	    names(X), "\n")		
	  cat("\n")
    options(warn=-0)
    }
	screened.parms <- names(X)[-1]
 
	if(p.sig) {
	cat("****  Evaluating and removing weak parameter contribution ", "\n")
	  flush.console(); Sys.sleep(0.01)
	  if(nrow(X) > 50000) {
	    X <- X[c(sample(which(X$treatment == 1), 25000), 
	             sample(which(X$treatment == 0), 25000)),] 
	  }  	  
	  X$treatment <- factor(X$treatment)
	  rf <- ranger(treatment ~ ., data = X, importance = "permutation") 
	  imp <- sort(importance(rf) / max(importance(rf)), decreasing = TRUE)
	      drop.imp <- names(imp)[which(imp < 0.001)]
		if(length(drop.imp) > 0)
		  screened.parms <- screened.parms[-which(screened.parms %in% drop.imp)]
		remove(rf)
	}
	remove(X)
  gc(verbose = FALSE, full = TRUE)

  #################################   
  #################################  
  if(!inherits(mdata, "list")){
    sub.models = FALSE
  	st <- system.time({
    cat("****   Model has", nrow(mdata), "observations", "\n") 
	cat("****     treatment n =", length(which(mdata$treatment == 1)), "\n") 
	cat("****     control n =", length(which(mdata$treatment == 0)), "\n")
	cat("****  started at", as.character(Sys.time()), "\n")	
    cat("***************************************************", "\n")	
    flush.console(); Sys.sleep(0.01)	

	X = mdata[,screened.parms]
    Y = mdata$y 
	W = as.numeric(mdata$treatment)
      cl <- makePSOCKcluster(15)   
        registerDoParallel(cl) 
          tau.forest <- causal_forest(
            X = X, Y = Y,
            W = W, num.trees = nboot,
            seed = 42, num.threads = nt)
	  stopCluster(cl)
    })	  
    #*************************************************
    cat("  Finished at", as.character(Sys.time()), "\n")
    cat("\n")   
	st <- st[3]/60
      if(st > 60) st <- round(st / 60, 4)	
	    cat("System run time", st, "\n")

  ################################# 		
  #################################  
  } else if(inherits(mdata, "list")){
	sub.models = TRUE
    cat("\n")
    cat("****     Full-model has", length(mdata)*mdl.size, "observations", "\n") 	
	cat("****     Subsample models have", mdl.size, "observations", "\n") 
    cat("****  started at", as.character(Sys.time()), "\n")	
      flush.console(); Sys.sleep(0.01)	  
	st <- system.time({
	  tau.forest <- list()
      tau.var <- vector()	  
	  for(d in 1:length(mdata)) {
  	    X = mdata[[d]][,screened.parms]
        Y = mdata[[d]]$y 
        W = as.numeric(mdata[[d]]$treatment)
		cat("****  Sub-model:", d, "of", length(mdata), "n = ", length(Y), "(", as.character(Sys.time()), ")", "\n")
          flush.console(); Sys.sleep(0.01)
        cl <- makeCluster(detectCores()-1, type = "SOCK") 
          registerDoParallel(cl) 
            tau.forest[[d]] <- causal_forest(
              X = X, Y = Y, 
              W = W, num.trees = nboot,
      		  ci.group.size = 2, seed = 42, 
			  num.threads = nt)
          stopCluster(cl)
		    tau.var <- append(tau.var, predict(object = tau.forest[[d]], 
			                  estimate.variance = TRUE)$variance.estimates)  
		  remove(X,Y,W)
        gc()
	  }
	  
	cat("**** Merging forests", "\n")
	cat("***************************************************", "\n")
      flush.console(); Sys.sleep(0.01)
      tau.forest <- merge_forests(tau.forest, compute.oob.predictions = FALSE)
    })
	st <- st[3]/60
      if(st > 60) st <- round(st / 60, 4)	
	    cat("System run time", st, "\n")
	cat("\n") 
      cat("  Finished at", as.character(Sys.time()), "\n")
        cat("\n")   
    }

  #################################
  #################################  
  # Write model results:
	ate <- average_treatment_effect(tau.forest, target.sample = "overlap", method = "TMLE")    
	  if(length(ate) == 2){
        cat("**** CATE for", j, ate[1], "\n")
        if(ate[1] > 1) ate[1] <- 1  
	    if(ate[1] < -1) ate[1] <- -1
	  }
	if(!exists("i")) i = 1  
	if(group.by == "pa") {
      cf.ate[[paste(j,i,sep="_")]] <- data.frame(NAME = j, effect.size = ate[1], std.error = ate[2])
    } else if(group.by == "elu"){
      cf.ate[[paste(j,i,sep="_")]] <- data.frame(elu = j, effect.size=ate[1], std.error = ate[2])
    }
	cat("\n")
    cat("**** Estimating Heterogeneous Treatment Effect Estimates for", j, "\n")
      flush.console(); Sys.sleep(0.01)
	st <- system.time({  
	if(!sub.models){
      preds <- predict(object = tau.forest, estimate.variance = TRUE)
	    preds <- preds[which(mdata$treatment == 1),]
          preds <- data.frame(cell = mdata$cell[which(mdata$treatment == 1)], preds)
	} else if(sub.models){
	  mdata <- do.call(rbind, mdata)
	    treat.idx <- which(mdata$treatment == 1) 
	      mdata <- mdata[treat.idx,]
	        preds <- predict(object = tau.forest, newdata = mdata[,screened.parms], 
		                     estimate.variance = TRUE, num.threads = nt/2)
	          preds <- data.frame(cell = mdata$cell, predictions = preds$predictions, 
			                      variance.estimates = tau.var[treat.idx])
	}
	nidx <- which(names(preds) %in% c("debiased.error", "excess.error"))
	  if(length(nidx) > 0) preds <- preds[,-nidx]
 	preds$unscaled.es <- preds$predictions 
	  pos.idx <- which(preds$predictions > 0)
	    if(length(pos.idx) > 0) 
          preds$predictions[pos.idx] <- preds$predictions[pos.idx] / max(preds$predictions[pos.idx])   
	  neg.idx <- which(preds$predictions < 0)
	    if(length(neg.idx) > 0) 
	      preds$predictions[neg.idx] <- preds$predictions[neg.idx] / min(preds$predictions[neg.idx]) * -1   	  
    didx <- which(duplicated(preds$cell)) 
      if(length(didx) > 0) { 
	  cat("**** Reconciling", length(didx), "duplicates", "\n")
        preds <- as.data.table(preds)
        preds <- as.data.frame(preds[,lapply(.SD, median), "cell"])
  	  }
	if("variance.estimates" %in% names(preds)){
	  preds$LCI <- preds$predictions - 1.96 * sqrt(preds$variance.estimates)
      preds$UCI <- preds$predictions + 1.96 * sqrt(preds$variance.estimates)
    }
	cat("**** Writing Heterogeneous Treatment Effect Estimates to raster", "\n")
	  if(nlyr(es.raster) == 5) {
	    es.raster[[1:5]][preds$cell] <- preds[,which(names(preds) %in% c("predictions", "variance.estimates", "unscaled.es", "LCI", "UCI"))]
	  } else if(nlyr(es.raster) == 2) {
        es.raster[[1:2]][preds$cell] <- preds[,which(names(preds) %in% c("predictions", "unscaled.es"))]
	  }
	  stars::write_stars(st_as_stars(es.raster), dsn = es.name,
                         type = "Float32", update = TRUE)
	})
	st <- st[3] / 60
      if(st > 60) st <- round(st / 60, 4)	
	    cat("System run time", st, "\n")
        if(write.db) {
		  cat("**** Writing effect sizes to database", "\n")
            preds <- data.frame(NAME = j, response = y, preds) 
          db.name <- file.path(pred.dir, paste(metric.type, mres, y, "effect_size.duckdb", sep="_"))
          con <- dbConnect(duckdb::duckdb(dbdir = db.name), read_only = FALSE)
            if(!DBI::dbExistsTable(con, "effect_sizes")) {
              duckdb::dbWriteTable(con, "effect_sizes", preds)
            } else {
              duckdb::dbAppendTable(con, "effect_sizes", preds) 
              duckdb::dbAppendTable(con, "effect_sizes", preds) 
            }
			dbDisconnect(con, shutdown = TRUE)
        }		
    remove(mdata, preds, ate, pos.idx, neg.idx, didx)  
  gc(verbose = FALSE, full = TRUE)   
} # End of model(s) for loop

cf.ate <- do.call(rbind, cf.ate)
  rownames(cf.ate) <- 1:nrow(cf.ate)
write.csv(cf.ate, file.path(pred.dir, paste(y, mres, "ATE.csv", sep="_")))


# plot results
clrs <- c("darkred","darkorange1","azure3","aquamarine3","darkgreen")
bks <- c( seq(-1,-0.005, 0.10), c(-0.01, 0.01), seq(0.10, 1, 0.1)) 
hte.mcls <- colorRampPalette(clrs)(21)      
plot(es.raster[[1]], maxcell=5000000, smooth = TRUE, type="interval",  
     breaks = bks, col = hte.mcls)
  plot(st_geometry(pa), add=TRUE)
    plot(st_geometry(pfp), add=TRUE)
      title(paste0(y, " effect size"))

#**************************************
#**************************************
