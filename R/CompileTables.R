# Summarize effect size and trend results into single spreadsheet(s)
root = "C:/evans/PFP/results"
cty.names = c("brazil", "bhutan", "canada", "colombia", "costa_rica", "peru")
table.type = c("interventions", "iucn", "pa", "pfp")[4]
cty = file.path(root, cty.names, "tables", table.type)  
  if(any(table.type %in% c("pa", "pfp"))) {
    cty <- cty[grep(paste(c("brazil", "colombia", "peru"), collapse="|"), cty)]
	cty.names = c("brazil", "colombia", "peru")
  }

#********************************************************
# Summarize effect sizes
es.names = c("metric", "class", "count", "ha", "pct")
mnames <- c("change_lt_0", "change_gt_0", "change_lt_neg0.7", "change_gt_0.7",
            "current_lt_0", "current_gt_0", "current_lt_neg0.7", "current_gt_0.7",
            "tsa_lt_0", "tsa_gt_0", "tsa_lt_neg0.7", "tsa_gt_0.7")

effect.sizes.lai <- data.frame(matrix(ncol = 13, nrow = length(cty.names)))
  names(effect.sizes.lai) <- c("country", mnames)
    effect.sizes.lai["country"] <- cty.names   
effect.sizes.fcov <- data.frame(matrix(ncol = 13, nrow = length(cty.names)))
  names(effect.sizes.fcov) <- c("country", mnames)
    effect.sizes.fcov["country"] <- cty.names   

for(d in cty) {
  setwd(d)
  cty.idx <- which(cty.names %in% unlist(strsplit(d, "/"))[5]) 
  es <- read.csv("effect_size_class_freq.csv")
    es <- split(es, f=es$metric)
	  es.lai <- es[grep("lai", names(es))]
	  es.fcov <- es[grep("fcov", names(es))]
    for(m in c("current", "change","tsa")) {
      # "metric_lt_0", "metric_gt_0", "metric_lt_neg07", "metric_gt_07" 
	  # Less than -0.70,  Greater than 0.7,  less than 0,  Greater than 0 		
	  
	  # LAI
	  dat <- es.lai[[grep(m, names(es.lai))]]  
        midx <- grep(m, names(effect.sizes.lai))
	  effect.sizes.lai[cty.idx,][midx] <- c(
        # pre_lt_0       
 	    paste0(format(round(sum(dat[1:10,]$ha),0), big.mark=",", scientific=FALSE, trim=TRUE) , " ha ",    
               paste0("(", paste0(round(sum(dat[1:10,]$pct) * 100, 4), "%"), ")")  ),  
	    # pre_gt_0
	    paste0(format(round(sum(dat[12:21,]$ha),0), big.mark=",", scientific=FALSE, trim=TRUE), " ha ", 
               paste0("(", paste0(round(sum(dat[12:21,]$pct) * 100, 4), "%"), ")")  ),
        # pre_lt_neg0.7  
	    paste0(format(round(sum(dat[1:3,]$ha),0), big.mark=",", scientific=FALSE, trim=TRUE), " ha ", 
               paste0("(", paste0(round(sum(dat[1:3,]$pct) * 100, 4), "%"), ")")  ),
		# pre_gt_0.7
 	    paste0(format(round(sum(dat[19:21,]$ha),0), big.mark=",", scientific=FALSE, trim=TRUE), " ha ", 
               paste0("(", paste0(round(sum(dat[19:21,]$pct) * 100, 4), "%"), ")")  ) )

      # FCOV
	  dat <- es.fcov[[grep(m, names(es.fcov))]]  
      midx <- grep(m, names(effect.sizes.fcov))	  
	  effect.sizes.fcov[cty.idx,][midx] <- c(
        # pre_lt_0       
 	    paste0(
        format(round(sum(dat[1:10,]$ha),0), big.mark=",", scientific=FALSE, trim=TRUE) , " ha ",    
        paste0("(", paste0(round(sum(dat[1:10,]$pct) * 100, 4), "%"), ")")  ),  
	    # pre_gt_0
	    paste0(	
        format(round(sum(dat[12:21,]$ha),0), big.mark=",", scientific=FALSE, trim=TRUE), " ha ", 
        paste0("(", paste0(round(sum(dat[12:21,]$pct) * 100, 4), "%"), ")")  ),
        # pre_lt_neg0.7  
	    paste0(	
        format(round(sum(dat[1:3,]$ha),0), big.mark=",", scientific=FALSE, trim=TRUE), " ha ", 
        paste0("(", paste0(round(sum(dat[1:3,]$pct) * 100, 4), "%"), ")")  ),
		# pre_gt_0.7
 	    paste0(	
        format(round(sum(dat[19:21,]$ha),0), big.mark=",", scientific=FALSE, trim=TRUE), " ha ", 
        paste0("(", paste0(round(sum(dat[19:21,]$pct) * 100, 4), "%"), ")")  ) )
  	}
}	
effect.sizes.lai <- data.frame(metric="LAI", effect.sizes.lai)
effect.sizes.fcov <- data.frame(metric="fCOV", effect.sizes.fcov)
effect.sizes <- rbind(effect.sizes.lai, effect.sizes.fcov)
write.csv(effect.sizes, file.path(root, "tables", paste0(table.type, "_EffectSizeSummaries.csv")))

#********************************************************
# Summarize trends	  
trend.names = c("period", "class", "count", "ha", "pct")
tnames <- c("pre_lt_0", "pre_gt_0", "pre_lt_neg0.7", "pre_gt_0.7",
            "post_lt_0", "post_gt_0", "post_lt_neg0.7", "post_gt_0.7")

trend.lai <- data.frame(matrix(ncol = 9, nrow = length(cty.names)))
  names(trend.lai) <- c("country", tnames)
    trend.lai["country"] <- cty.names 
	
trend.fcov <- data.frame(matrix(ncol = 9, nrow = length(cty.names)))
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
	  trend.lai[cty.idx,][midx] <- c(
        # pre_lt_0       
 	    paste0(
        format(round(sum(dat[1:10,]$ha),0), big.mark=",", scientific=FALSE, trim=TRUE) , " ha ",    
        paste0("(", paste0(round(sum(dat[1:10,]$pct) * 100, 4), "%"), ")")  ),  
	    # pre_gt_0
	    paste0(	
        format(round(sum(dat[12:21,]$ha),0), big.mark=",", scientific=FALSE, trim=TRUE), " ha ", 
        paste0("(", paste0(round(sum(dat[12:21,]$pct) * 100, 4), "%"), ")")  ),
        # pre_lt_neg0.7  
	    paste0(	
        format(round(sum(dat[1:3,]$ha),0), big.mark=",", scientific=FALSE, trim=TRUE), " ha ", 
        paste0("(", paste0(round(sum(dat[1:3,]$pct) * 100, 4), "%"), ")")  ),
		# pre_gt_0.7
 	    paste0(	
        format(round(sum(dat[19:21,]$ha),0), big.mark=",", scientific=FALSE, trim=TRUE), " ha ", 
        paste0("(", paste0(round(sum(dat[19:21,]$pct) * 100, 4), "%"), ")")  ) )

	  # fCOV
	  dat <- tfcov[[grep(m, names(tfcov))]]  
      midx <- grep(m, names(trend.fcov))	  
	  trend.fcov[cty.idx,][midx] <- c(
        # pre_lt_0       
 	    paste0(
        format(round(sum(dat[1:10,]$ha),0), big.mark=",", scientific=FALSE, trim=TRUE) , " ha ",    
        paste0("(", paste0(round(sum(dat[1:10,]$pct) * 100, 4), "%"), ")")  ),  
	    # pre_gt_0
	    paste0(	
        format(round(sum(dat[12:21,]$ha),0), big.mark=",", scientific=FALSE, trim=TRUE), " ha ", 
        paste0("(", paste0(round(sum(dat[12:21,]$pct) * 100, 4), "%"), ")")  ),
        # pre_lt_neg0.7  
	    paste0(	
        format(round(sum(dat[1:3,]$ha),0), big.mark=",", scientific=FALSE, trim=TRUE), " ha ", 
        paste0("(", paste0(round(sum(dat[1:3,]$pct) * 100, 4), "%"), ")")  ),
		# pre_gt_0.7
 	    paste0(	
        format(round(sum(dat[19:21,]$ha),0), big.mark=",", scientific=FALSE, trim=TRUE), " ha ", 
        paste0("(", paste0(round(sum(dat[19:21,]$pct) * 100, 4), "%"), ")")  ) )
  }
}	
trend.lai <- data.frame(metric="LAI", trend.lai)
trend.fcov <- data.frame(metric="fCOV", trend.fcov)
trends <- rbind(trend.lai, trend.fcov)
write.csv(trends, file.path(root, "tables",paste0(table.type, "_TrendsSummaries.csv")))
