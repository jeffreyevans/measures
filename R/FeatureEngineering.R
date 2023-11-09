# Feature engineering to produce design matrix covariates
#   that represent characteristics of the timeseries
# https://cran.r-project.org/web/packages/tsfeatures/vignettes/tsfeatures.html

#x <- as.numeric(intervention[1,][8:ncol(intervention)])

#*************************************************
# Feature Engineering function
feature.engineering <- function(tsdata, ID=NULL, year.month = c(2014, 1), 
                                ts.freq = 36, trace = TRUE) {
  suppressMessages(lapply(c("dplyr", "lubridate", "tsfeatures"), 
	  	           require, character.only = TRUE))
  if(is.null(ID)) ID <- 1:nrow(tsdata) 
  process.ts <- function(j, start = year.month, f=ts.freq) {
    ts.features <- function(x) {
      tsf <- c(entropy(x), stability(x), lumpiness(x),
        hurst(x), max_level_shift(x), max_var_shift(x),
        max_kl_shift(x), crossing_points(x),
        flat_spots(x), stl_features(x)[c(3,4,6,7,8)],
        acf_features(x), pacf_features(x), arch_stat(x),
        heterogeneity(x), firstmin_ac(x))
      return(tsf)
    } 
    j <- ts(j, start=start, frequency=f)
        ts.metrics <- ts.features(j)  
      names(ts.metrics)[34] <- "facf"  
    return(ts.metrics)
  }  
  heterogeneity <- function(x) {
      x.whitened <- na.contiguous(ar(x)$resid)
      x.archtest <- arch_stat(as.matrix(x.whitened))
      LBstat <- sum(acf(x.whitened^2, lag.max = 12L, plot = FALSE)$acf[-1L]^2)
      garch.fit <- suppressWarnings(tseries::garch(x.whitened, trace = FALSE))
      garch.fit.std <- residuals(garch.fit)
      x.garch.archtest <- arch_stat(garch.fit.std)
      LBstat2 <- NA
      try(LBstat2 <- sum(acf(na.contiguous(garch.fit.std^2), lag.max = 12L, 
          plot = FALSE)$acf[-1L]^2), silent = TRUE)
      output <- c(arch_acf = LBstat, garch_acf = LBstat2, arch_r2 = unname(x.archtest), 
          garch_r2 = unname(x.garch.archtest))
      return(output)
  }   
  tfe <- setNames(vector("list", nrow(tsdata)), ID)
    cf.idx <- ID
    for(i in 1:length(tfe)) {
	  if(trace)
	    cat("processing", i, "of", length(tfe), "\n")
	  v <- cf.idx[i]
      tryCatch({
        fe <- process.ts(na.omit(as.numeric(tsdata[i,])))
          }, error=function(e){ fe <- rep(NA,34) })
      tfe[[as.character(v)]] <- fe
    }	
      tfe <- do.call(rbind, tfe)
    tfe <- data.frame(ID=ID, tfe)
  return(tfe)	
} 
