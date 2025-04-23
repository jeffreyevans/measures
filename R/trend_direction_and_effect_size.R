library(terra)
setwd("C:/evans/PFP/results/costa_rica")

lai <- rast("lai_effect_sizes_classified.tif")
lai.trend <- mask(rast("lai_tau_change.tif"),lai[[1]])

d <- na.omit(data.frame(c(lai,lai.trend)))
  ft <- lapply(c("current", "gain", "loss", "tsa"), \(i) { 
    f <- data.frame(metric = i, table(d[,"lai_change_type"], d[,i]))
      names(f) <- c("metric",  "trend", "effect_size", "count")
  	  return(f)
  })
f <- do.call("rbind", ft)

write.csv(f, file.path(getwd(), "tables", "effect_size_trend_frequencies.csv"), row.names = FALSE)
