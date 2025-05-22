library(terra)

root = "C:/evans/PFP/results"
setwd(root)
country <- c("brazil", "bhutan", "canada", "colombia", "costa_rica", "peru")

ss <- data.frame(matrix(ncol = 4, nrow =length(country) ))
  names(ss) <- c("country", "treatment_area", "control_area", "sample_pct")
    ss[,1] <- country
	  ss[,1][which(ss[,1] %in% c("brazil", "canada"))] <- c("ARPA", "Great Bear")

  for(i in 1:length(country)) {
    treat <- rast( file.path(root,  country[i], "lai_treatment_responses.tif") ) 
      ss[i,][2] <- expanse(treat[[3]], byValue=FALSE, transform=FALSE, usenames=TRUE, unit="ha")$area
    ctl <- rast( file.path(root,  country[i], "lai_control_responses.tif") )   
      ss[i,][3] <- expanse(ctl[[3]], byValue=FALSE, transform=FALSE, usenames=TRUE, unit="ha")$area
  }

ss[,"sample_pct"] <- ss$control_area / ss$treatment_area 

