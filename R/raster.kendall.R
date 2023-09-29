#' @param x             A multiband terra SpatRaster object with at least 5 layers
#' @param intercept     (FALSE/TRUE) return a raster with the pixel 
#'                       wise intercept values 
#' @param p.value       (FALSE/TRUE) return a raster with the pixel 
#'                       wise p.values
#' @param z.value       (FALSE/TRUE) return a raster with the pixel 
#'                       wise z.values
#' @param confidence    (FALSE/TRUE) return a raster with the pixel 
#'                       wise 95 pct confidence levels
#' @param tau           (FALSE/TRUE) return a raster with the pixel wise 
#'                       tau correlation values
#' @param ...           Additional arguments passed to the raster 

kendall <- function(y, tau = TRUE, intercept = TRUE, p.value = TRUE, 
                    z.value = TRUE, confidence = TRUE, 
					prewhiten = FALSE, na.rm, ...) {
    if(length(y[!is.na(y)]) < 8) 
      stop("The Kendall Tau needs at least 8 observations")
    pass.sum <- 0
	out.names <- 
	c("slope", "tau", "intercept", "p-value", "z-value", "limits.LCL", "limits.UCL")[
	  which(c(TRUE, tau, intercept, p.value, z.value,  rep(confidence,2)))]
	if(prewhiten) {
	  confidence = FALSE
	  intercept = FALSE 
    }
      if( p.value ) pass.sum = pass.sum + 1
	    if( z.value ) pass.sum = pass.sum + 1
	    if( tau ) pass.sum = pass.sum + 1
	  if( confidence ) pass.sum = pass.sum + 2
	if( intercept ) pass.sum = pass.sum + 1
      fit.results <- c(rep(NA,pass.sum + 1))
    if(!prewhiten) {
      if(!any(which(utils::installed.packages()[,1] %in% "EnvStats")))
        stop("please install EnvStats package before running this function")	
      fit <- EnvStats::kendallTrendTest(y ~ 1)
      fit.results <- fit$estimate[2]
        if(tau == TRUE) { fit.results <- c(fit.results, fit$estimate[1]) }
            if(intercept == TRUE) { fit.results <- c(fit.results, fit$estimate[3]) }  
              if(p.value == TRUE) { fit.results <- c(fit.results, fit$p.value) } 
                if(z.value == TRUE) { fit.results <- c(fit.results, fit$statistic) }
	        if(confidence == TRUE) { 
          ci <- unlist(fit$interval["limits"])
            if( length(ci) == 2) { 
              fit.results <- c(fit.results, ci)
              } else {
                fit.results <- c(fit.results, c(NA,NA))
              }			  
         }
    } else {
    # kendall autocorrelation correction (pre-whitening)
    x = y
	  z = NULL
        pval = NULL
        S = 0
      var.S = NULL
    Tau = NULL
      if (any(is.finite(x) == FALSE)) {
        x <- x[-c(which(is.finite(x) == FALSE))]
      }
      n <- length(x)
    V <- rep(NA, n * (n - 1)/2)
    k = 0
      for (i in 1:(n - 1)) {
        for (j in (i + 1):n) {
          k = k + 1
          V[k] = (x[j] - x[i])/(j - i)
        }
      }
    slp <- stats::median(V, na.rm = TRUE)
      t1 = 1:length(x)
        xt <- (x[1:n]) - ((slp) * (t1))
    ro <- stats::acf(xt, lag.max = 1, plot = FALSE)$acf[-1]
      a = 1:(length(xt) - 1)
        b = 2:(length(xt))
          xp <- (xt[b] - (xt[a] * ro))
          l <- length(xp)
        q = 1:l
      y <- (xp[1:l] + ((slp) * (q)))
    n1 <- length(y)
      for (i in 1:(n1 - 1)) {
        for (j in (i + 1):n1) {
          S = S + sign(y[j] - y[i])
        }
      }
    var.S = n1 * (n1 - 1) * (2 * n1 + 5) * (1/18)
      if (length(unique(y)) < n1) {
        aux <- unique(y)
          for (i in 1:length(aux)) {
            tie <- length(which(y == aux[i]))
              if (tie > 1) {
                var.S = var.S - tie * (tie - 1) * (2 * tie + 5) * (1/18)
              }
          }
      }
    if (S == 0) { z = 0 }
    if (S > 0) {
      z = (S - 1)/sqrt(var.S)
    } else {
      z = (S + 1)/sqrt(var.S)
    }
    pval = 2 * stats::pnorm(-abs(z))
      Tau = S/(0.5 * n1 * (n1 - 1))
        W <- rep(NA, n1 * (n1 - 1)/2)
          m = 0
    for (i in 1:(n1 - 1)) {
      for (j in (i + 1):n1) {
        m = m + 1
        W[m] = (y[j] - y[i])/(j - i)
      }
    }
    slp1 <- stats::median(W, na.rm = TRUE)
	  fit.results <- slp1 
	  if(tau == TRUE) { fit.results <- c(fit.results, Tau) }
        if(p.value == TRUE) { fit.results <- c(fit.results, pval) } 
          if(z.value == TRUE) { fit.results <- c(fit.results, z) }
	  fit.results <- as.numeric(fit.results) 
	    
    }
	names(fit.results) <- out.names
  return(fit.results)
}	
 
raster.kendall <- function(x, intercept = TRUE, p.value = TRUE, z.value = TRUE,   
                           confidence = TRUE, tau = TRUE, ...) {
  if(!any(which(utils::installed.packages()[,1] %in% "EnvStats")))
    stop("please install EnvStats package before running this function")
   if (!inherits(x, "SpatRaster")) 
	  stop(deparse(substitute(x)), " must be a terra SpatRaster object")
  if(confidence) {confidence = c(TRUE,TRUE)} else {confidence = c(FALSE,FALSE)}
    n <- c("intercept", "p.value", "z.value", "LCI", "UCI", "tau")	
	n <- n[which(c(intercept, p.value, z.value,confidence, tau))]	
    if( terra::nlyr(x) < 5) 
	  stop("Too few layers (n < 5) to calculate a trend")
  trend.slope <- function(y, tau.pass = tau, p.value.pass = p.value,  
                          confidence.pass = confidence[1], z.value.pass = z.value,
                          intercept.pass = intercept) {
    fit <- suppressWarnings( EnvStats::kendallTrendTest(y ~ 1) )
      fit.results <- fit$estimate[2]
        if(p.value.pass == TRUE) { fit.results <- c(fit.results, fit$p.value) } 
		  if(z.value.pass == TRUE) { fit.results <- c(fit.results, fit$statistic) } 
  	        if(confidence.pass == TRUE) { 
		      ci <- unlist(fit$interval["limits"])
		        if( length(ci) == 2) { 
		          fit.results <- c(fit.results, ci)
                } else {
                  fit.results <- c(fit.results, c(NA,NA))
                }			  
		    }
        if(intercept.pass == TRUE) { fit.results <- c(fit.results, fit$estimate[3]) }  
		  if(tau.pass == TRUE) { fit.results <- c(fit.results, fit$estimate[1]) }  
    return( fit.results )
  }
  k <- terra::app(x, fun=trend.slope, ...)
    names(k) <- c("slope", n)
  return( k )
}
