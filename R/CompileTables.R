# Summarize effect size and trend results into single spreadsheet(s)
root = "C:/evans/PFP/results"
cty.names = c("brazil", "bhutan", "canada", "colombia", "costa_rica", "peru")

table.type = c("interventions", "iucn", "pa", "pfp")[4]
cty = file.path(root, cty.names, "tables", table.type)  
  if(any(table.type %in% c("pa", "pfp"))) {
    cty <- cty[grep(paste(c("brazil", "colombia", "peru"), collapse="|"), cty)]
	cty.names = c("brazil", "colombia", "peru")
  }

run.es = c(TRUE, FALSE)[1]
run.trends = c(TRUE, FALSE)[1]

#********************************************************
# Summarize effect sizes
if(run.es) {

es.names = c("metric", "class", "count", "ha", "pct")
mnames <- c("change_lt_neg0.7", "change_lt_0", "change_zero", "change_small", "change_gt_0",  "change_gt_0.7",
            "current_lt_neg0.7", "current_lt_0", "current_zero", "current_small", "current_gt_0",  "current_gt_0.7",
            "tsa_lt_neg0.7", "tsa_lt_0", "tsa_zero", "tsa_small", "tsa_gt_0",  "tsa_gt_0.7")

effect.sizes.lai <- data.frame(matrix(ncol = length(mnames)+1, nrow = length(cty.names)))
  names(effect.sizes.lai) <- c("country", mnames)
    effect.sizes.lai["country"] <- cty.names   
effect.sizes.fcov <- data.frame(matrix(ncol = length(mnames)+1, nrow = length(cty.names)))
  names(effect.sizes.fcov) <- c("country", mnames)
    effect.sizes.fcov["country"] <- cty.names   

for(d in cty) {
  setwd(d)
  cty.idx <- which(cty.names %in% unlist(strsplit(d, "/"))[5]) 
  es <- read.csv("effect_size_class_freq.csv")
    es <- split(es, f=es$metric)
	  es.lai <- es[grep("lai", names(es))]
	  es.fcov <- es[grep("fcov", names(es))]
    for(m in c("change", "current", "tsa")) {
	  # LAI
	  dat <- es.lai[[grep(m, names(es.lai))]]  
        midx <- grep(m, names(effect.sizes.lai))
          idx0 <- which(dat$class %in% c("-0.01 to 0.01"))	
	      idx1 <- which(dat$class %in% c("-1 to -0.9", "-0.9 to -0.8", "-0.8 to -0.7")) 	
	      idx2 <- grep("-", dat$class)
          idx3 <- which(dat$class %in% c("-0.1 to -0.01", "-0.01 to 0.01", "0.01 to 0.1"))
          idx4 <- which(!grepl("-", dat$class))	
	      idx5 <- which(dat$class %in% c("0.7 to 0.8", "0.8 to 0.9", "0.9 to 1"))
	  effect.sizes.lai[cty.idx,][midx] <- c(
        # pre_lt_neg0.7
	    paste0(format(round(sum(dat[idx1,]$ha),0), big.mark=",", scientific=FALSE, trim=TRUE), " ha ", 
               paste0("(", paste0(round(sum(dat[idx1,]$pct) * 100, 4), "%"), ")")  ),	  
        # pre_lt_0
 	    paste0(format(round(sum(dat[idx2,]$ha),0), big.mark=",", scientific=FALSE, trim=TRUE) , " ha ",    
               paste0("(", paste0(round(sum(dat[idx2,]$pct) * 100, 4), "%"), ")")  ),  
		# zero
 	    paste0(format(round(sum(dat[idx0,]$ha),0), big.mark=",", scientific=FALSE, trim=TRUE) , " ha ",    
               paste0("(", paste0(round(sum(dat[idx0,]$pct) * 100, 4), "%"), ")")  ),  		
        # small values (-0.1 - 0.1)
 	    paste0(format(round(sum(dat[idx3,]$ha),0), big.mark=",", scientific=FALSE, trim=TRUE) , " ha ",    
               paste0("(", paste0(round(sum(dat[idx3,]$pct) * 100, 4), "%"), ")")  ),  		
	    # pre_gt_0
	    paste0(format(round(sum(dat[idx4,]$ha),0), big.mark=",", scientific=FALSE, trim=TRUE), " ha ", 
               paste0("(", paste0(round(sum(dat[idx4,]$pct) * 100, 4), "%"), ")")  ),
		# pre_gt_0.7
 	    paste0(format(round(sum(dat[idx5,]$ha),0), big.mark=",", scientific=FALSE, trim=TRUE), " ha ", 
               paste0("(", paste0(round(sum(dat[idx5,]$pct) * 100, 4), "%"), ")")  ) )

      # FCOV
	  dat <- es.fcov[[grep(m, names(es.fcov))]]  
        midx <- grep(m, names(effect.sizes.fcov))
          idx0 <- which(dat$class %in% c("-0.01 to 0.01"))		
	      idx1 <- which(dat$class %in% c("-1 to -0.9", "-0.9 to -0.8", "-0.8 to -0.7")) 	
	      idx2 <- grep("-", dat$class)
          idx3 <- which(dat$class %in% c("-0.1 to -0.01", "-0.01 to 0.01", "0.01 to 0.1"))
          idx4 <- which(!grepl("-", dat$class))	
	      idx5 <- which(dat$class %in% c("0.7 to 0.8", "0.8 to 0.9", "0.9 to 1"))  
	  effect.sizes.fcov[cty.idx,][midx] <- c(
		# pre_lt_neg0.7  		
	    paste0(	
        format(round(sum(dat[idx1,]$ha),0), big.mark=",", scientific=FALSE, trim=TRUE), " ha ", 
        paste0("(", paste0(round(sum(dat[idx1,]$pct) * 100, 4), "%"), ")")  ),
        # pre_lt_0 	
 	    paste0(
        format(round(sum(dat[idx2,]$ha),0), big.mark=",", scientific=FALSE, trim=TRUE) , " ha ",    
        paste0("(", paste0(round(sum(dat[idx2,]$pct) * 100, 4), "%"), ")")  ),          
		# zero
 	    paste0(format(round(sum(dat[idx0,]$ha),0), big.mark=",", scientific=FALSE, trim=TRUE) , " ha ",    
               paste0("(", paste0(round(sum(dat[idx0,]$pct) * 100, 4), "%"), ")")  ),  			    				
		# small values (-0.1 - 0.1)
 	    paste0(format(round(sum(dat[idx3,]$ha),0), big.mark=",", scientific=FALSE, trim=TRUE) , " ha ",    
               paste0("(", paste0(round(sum(dat[idx3,]$pct) * 100, 4), "%"), ")")  ),  			    
		
		# pre_gt_0
	    paste0(	
        format(round(sum(dat[idx4,]$ha),0), big.mark=",", scientific=FALSE, trim=TRUE), " ha ", 
        paste0("(", paste0(round(sum(dat[idx4,]$pct) * 100, 4), "%"), ")")  ),
		# pre_gt_0.7
 	    paste0(	
        format(round(sum(dat[idx5,]$ha),0), big.mark=",", scientific=FALSE, trim=TRUE), " ha ", 
        paste0("(", paste0(round(sum(dat[idx5,]$pct) * 100, 4), "%"), ")")  ) )
  	}
}	
  effect.sizes.lai <- data.frame(metric="LAI", effect.sizes.lai)
  effect.sizes.fcov <- data.frame(metric="fCOV", effect.sizes.fcov)
  effect.sizes <- rbind(effect.sizes.lai, effect.sizes.fcov)
  write.csv(effect.sizes, file.path(root, "tables", paste0(table.type, "_EffectSizeSummaries.csv")))
}

#********************************************************
# Summarize trends	 
if(run.trends) {
 
trend.names = c("period", "class", "count", "ha", "pct")
tnames <- c("pre_lt_neg0.7", "pre_lt_0", "pre_zero", "pre_small", "pre_gt_0",  "pre_gt_0.7",
            "post_lt_neg0.7", "post_lt_0", "post_zero", "post_small", "post_gt_0",  "post_gt_0.7")

trend.lai <- data.frame(matrix(ncol = length(tnames)+1, nrow = length(cty.names)))
  names(trend.lai) <- c("country", tnames)
    trend.lai["country"] <- cty.names 	
trend.fcov <- data.frame(matrix(ncol = length(tnames)+1, nrow = length(cty.names)))
  names(trend.fcov) <- c("country", tnames)
    trend.fcov["country"] <- cty.names   

for(d in cty) {
  setwd(d)
  cty.idx <- which(cty.names %in% unlist(strsplit(d, "/"))[5]) 
  trend <- read.csv("trend_class_freq.csv")
    trend <- split(trend, f=trend$period)
  	  tlai <- trend[grep("lai", names(trend))] 
  	  tfcov <- trend[grep("fcov", names(trend))] 
    for(m in c("pre", "post")) {
	  # LAI
	  dat <- tlai[[grep(m, names(tlai))]]  
      midx <- grep(m, names(trend.lai))	  
          idx0 <- which(dat$class %in% c("-0.01 to 0.01"))	
	      idx1 <- which(dat$class %in% c("-1 to -0.9", "-0.9 to -0.8", "-0.8 to -0.7")) 	
	      idx2 <- grep("-", dat$class)
          idx3 <- which(dat$class %in% c("-0.1 to -0.01", "-0.01 to 0.01", "0.01 to 0.1"))
          idx4 <- which(!grepl("-", dat$class))	
	      idx5 <- which(dat$class %in% c("0.7 to 0.8", "0.8 to 0.9", "0.9 to 1"))
	  trend.lai[cty.idx,][midx] <- c(
        # pre_lt_neg0.7
	    paste0(format(round(sum(dat[idx1,]$ha),0), big.mark=",", scientific=FALSE, trim=TRUE), " ha ", 
               paste0("(", paste0(round(sum(dat[idx1,]$pct) * 100, 4), "%"), ")")  ),	  
        # pre_lt_0
 	    paste0(format(round(sum(dat[idx2,]$ha),0), big.mark=",", scientific=FALSE, trim=TRUE) , " ha ",    
               paste0("(", paste0(round(sum(dat[idx2,]$pct) * 100, 4), "%"), ")")  ),  
		# zero
 	    paste0(format(round(sum(dat[idx0,]$ha),0), big.mark=",", scientific=FALSE, trim=TRUE) , " ha ",    
               paste0("(", paste0(round(sum(dat[idx0,]$pct) * 100, 4), "%"), ")")  ),  		
        # small values (-0.1 - 0.1)
 	    paste0(format(round(sum(dat[idx3,]$ha),0), big.mark=",", scientific=FALSE, trim=TRUE) , " ha ",    
               paste0("(", paste0(round(sum(dat[idx3,]$pct) * 100, 4), "%"), ")")  ),  		
	    # pre_gt_0
	    paste0(format(round(sum(dat[idx4,]$ha),0), big.mark=",", scientific=FALSE, trim=TRUE), " ha ", 
               paste0("(", paste0(round(sum(dat[idx4,]$pct) * 100, 4), "%"), ")")  ),
		# pre_gt_0.7
 	    paste0(format(round(sum(dat[idx5,]$ha),0), big.mark=",", scientific=FALSE, trim=TRUE), " ha ", 
               paste0("(", paste0(round(sum(dat[idx5,]$pct) * 100, 4), "%"), ")")  ) )


	  # fCOV
	  dat <- tfcov[[grep(m, names(tfcov))]]  
      midx <- grep(m, names(trend.fcov))	  
          idx0 <- which(dat$class %in% c("-0.01 to 0.01"))	
	      idx1 <- which(dat$class %in% c("-1 to -0.9", "-0.9 to -0.8", "-0.8 to -0.7")) 	
	      idx2 <- grep("-", dat$class)
          idx3 <- which(dat$class %in% c("-0.1 to -0.01", "-0.01 to 0.01", "0.01 to 0.1"))
          idx4 <- which(!grepl("-", dat$class))	
	      idx5 <- which(dat$class %in% c("0.7 to 0.8", "0.8 to 0.9", "0.9 to 1"))
	  trend.fcov[cty.idx,][midx] <- c(
        # pre_lt_neg0.7
	    paste0(format(round(sum(dat[idx1,]$ha),0), big.mark=",", scientific=FALSE, trim=TRUE), " ha ", 
               paste0("(", paste0(round(sum(dat[idx1,]$pct) * 100, 4), "%"), ")")  ),	  
        # pre_lt_0
 	    paste0(format(round(sum(dat[idx2,]$ha),0), big.mark=",", scientific=FALSE, trim=TRUE) , " ha ",    
               paste0("(", paste0(round(sum(dat[idx2,]$pct) * 100, 4), "%"), ")")  ),  
		# zero
 	    paste0(format(round(sum(dat[idx0,]$ha),0), big.mark=",", scientific=FALSE, trim=TRUE) , " ha ",    
               paste0("(", paste0(round(sum(dat[idx0,]$pct) * 100, 4), "%"), ")")  ),  		
        # small values (-0.1 - 0.1)
 	    paste0(format(round(sum(dat[idx3,]$ha),0), big.mark=",", scientific=FALSE, trim=TRUE) , " ha ",    
               paste0("(", paste0(round(sum(dat[idx3,]$pct) * 100, 4), "%"), ")")  ),  		
	    # pre_gt_0
	    paste0(format(round(sum(dat[idx4,]$ha),0), big.mark=",", scientific=FALSE, trim=TRUE), " ha ", 
               paste0("(", paste0(round(sum(dat[idx4,]$pct) * 100, 4), "%"), ")")  ),
		# pre_gt_0.7
 	    paste0(format(round(sum(dat[idx5,]$ha),0), big.mark=",", scientific=FALSE, trim=TRUE), " ha ", 
               paste0("(", paste0(round(sum(dat[idx5,]$pct) * 100, 4), "%"), ")")  ) )
  }
}	
  trend.lai <- data.frame(metric="LAI", trend.lai)
  trend.fcov <- data.frame(metric="fCOV", trend.fcov)
  trends <- rbind(trend.lai, trend.fcov)
  write.csv(trends, file.path(root, "tables",paste0(table.type, "_TrendsSummaries.csv")))
}
