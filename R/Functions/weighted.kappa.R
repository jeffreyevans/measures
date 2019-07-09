#' @title Weighted Kappa
#' @description Derives the Weighted Kappa with optional confidence intervals
#' 
#' @param x             A matrix, table or vector (of y provided) 
#' @param y             (NULL) optional comparison vector equal to x
#' @param weights       A weights matrix (same dim as x) or "Unweighted", "Equal-Spacing", "Fleiss-Cohen" 
#' @param conf.level    A confidence level for the statistic, default is none
#'
#' @return A vector with the kappa statistic and optional confidence intervals. 
#'
#' @note ... 
#'
#' @author Jeffrey S. Evans  <jeffrey_evans@@tnc.org> 
#'
#' @references Cohen, J. (1960) A coefficient of agreement for nominal scales. Educational and 
#'               Psychological Measurement, 20, 37-46.
#' @references Everitt, B.S. (1968), Moments of statistics kappa and weighted kappa. 
#'               The British Journal of Mathematical and Statistical Psychology, 21, 97-103.
#' @references Fleiss, J.L., Cohen, J., and Everitt, B.S. (1969), Large sample standard 
#'               errors of kappa and weighted kappa. Psychological Bulletin, 72, 332-327.
#' 
#' @examples
#' n <- c("<10%", "11-20%", "21-30%", "31-40%", "41-50%", ">50%")
#' m <- matrix(c(5,8,1,2,4,2, 3,5,3,5,5,0, 1,2,6,11,2,1,
#'               0,1,5,4,3,3, 0,0,1,2,5,2, 0,0,1,2,1,4), nrow=6, byrow=TRUE,
#'               dimnames = list(r1 = n, r2 = n) )
#' weighted.kappa(m, weight="Equal-Spacing")
#' 
#' # supply an explicit weight matrix
#' ncol(m)
#' (wm <- outer(1:ncol(m), 1:ncol(m), function(x, y) {
#'         1 - ((abs(x-y)) / (ncol(m)-1)) } ))
#' weighted.kappa(m, weight=wm, conf.level=0.95)
#' 
#' # Fleiss weights the similarities
#' fleiss <- matrix(c(106,10,4,22,28,10,2,12,6), 
#'                  ncol=3, byrow=TRUE)
#' w <- matrix(c(1.0000, 0.0000, 0.4444,0.0000, 1.0000, 0.6666, 
#'             0.4444, 0.6666, 1.0000), ncol=3)
#' weighted.kappa(x = fleiss, weights = w, conf.level = 0.95)
#'
#' @ export weighted.kappa
weighted.kappa <- function (x, y = NULL, weights = c("Unweighted", "Equal-Spacing", 
                            "Fleiss-Cohen"), conf.level = NA, ...) {
    if (is.character(weights)) weights <- match.arg(weights)
    if (!is.null(y)) {
      if (!identical(weights, "Unweighted")) 
        stop("Vector interface for weighted Kappa is not supported. Provide confusion matrix.")
      x <- factor(x)
        y <- factor(y)
          lvl <- unique(c(levels(x), levels(y)))
          x <- factor(x, levels = lvl)
        y <- factor(y, levels = lvl)
      x <- table(x, y, ...)
    } else {
      d <- dim(x)
      if (d[1L] != d[2L]) 
        stop("x must be square matrix if provided as confusion matrix")
    }
    d <- diag(x)
      n <- sum(x)
        nc <- ncol(x)
        colFreqs <- colSums(x)/n
      rowFreqs <- rowSums(x)/n
    kappa <- function(po, pc) (po - pc)/(1 - pc)
    std <- function(po, pc, W = 1) sqrt(sum(W * W * po * (1 - po)) / crossprod(1 - pc)/n)
      po <- sum(d)/n
        pc <- crossprod(colFreqs, rowFreqs)
          k <- kappa(po, pc)
            s <- std(po, pc)
    if (is.matrix(weights)) { 
	  W <- weights 
    } else if(weights == "Equal-Spacing") { 
      W <- 1 - abs(outer(1:nc, 1:nc, "-"))/(nc - 1)
    } else { 
	  W <- 1 - (abs(outer(1:nc, 1:nc, "-"))/(nc - 1))^2
	}
    pow <- sum(W * x)/n
      pcw <- sum(W * colFreqs %o% rowFreqs)
        kw <- kappa(pow, pcw)
          sw <- std(x/n, 1 - pcw, W)
    if (is.na(conf.level)) {
      if (identical(weights, "Unweighted")) { 
        wk <- as.vector(k)
      } else { 
	    wk <- as.vector(kw)
	  }	 
    } else {
      if (identical(weights, "Unweighted")) {
        ci <- as.vector(k) + c(1, -1) * qnorm((1 - conf.level)/2) * as.vector(s)
        wk <- c(kappa = k, lwr.ci = ci[1], upr.ci = ci[2])
    } else {
      ci <- kw + c(1, -1) * qnorm((1 - conf.level)/2) * as.vector(sw)
      wk <- c(kappa = kw, lwr.ci = ci[1], upr.ci = ci[2])
    }
  }
  return(wk)
}
