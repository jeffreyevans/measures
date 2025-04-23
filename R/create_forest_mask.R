# https://www.sciencedirect.com/science/article/abs/pii/S0034425712002088
# 
#
#
suppressMessages(
  lapply(c("sf", "spatialEco", "terra"), 
		   require, character.only = TRUE))
terraOptions(tempdir = "C:/temp")

mres = c("250m", "300m")[1]                       # resoultion of trend data
metric.type <- c("lai", "fcov")[1]                # which trend metric to use
  if(mres == "250m") metric.type <- "lai"         # If res 250 the only option is lai
  
idx = 2  # indicates country and projection 
  country <- c("bhutan", "costa_rica", "canada", "colombia", "peru", "brazil")[idx]
root = file.path("C:/evans/PFP", country)
setwd(file.path(root, "data"))
  results.dir = file.path(dirname(root), "results", country)
  dat.dir = file.path(root, "data")

#*************************************************
# Set intervention date
if(country == "brazil") {
  intervention.date = as.Date("2014-01-10")  
} else if(country == "bhutan") {
  intervention.date = as.Date("2017-01-10")  
} else if(country == "canada") {
  intervention.date = as.Date("2006-01-10")
} else if(country == "colombia") {
  intervention.date = as.Date("2020-01-10")  	
} else if(country == "costa_rica") {
  intervention.date = as.Date("2010-01-10")  
} else if(country == "peru") {
  intervention.date = as.Date("2019-01-10")  
} else {
  warning("Not a valid option")
}  
if(mres == "300m" & intervention.date < as.Date("2015-01-10")) {
  intervention.date = as.Date("2015-01-10")
}

#*************************************************
# Parameters 
forest.threshold = c(2, 2.5, 1.5, 2.5, 2.5, 2.5)[idx]     # LAI forest threshold
forest.pct = 0.10                                         # fraction in the timeseries that lai >= forest.threshold 

#*************************************************
# Read metric, calculate max LAI and mask to forest
m <- rast(paste0("forest_", mres, ".tif"))
  metric <- rast(paste0(toupper(metric.type), "_", mres, "_trend.tif"))
    dates <- as.Date(names(metric))

pre.idx <- which(dates < intervention.date)   
post.idx <- which(dates >= intervention.date)   

full.fpct <- app(metric, fun=\(j) { mean(j >= forest.threshold) } ) 
  full.forest <- ifel(full.fpct >= forest.pct, 1, 0)
    freq(full.forest)
      full.forest[full.forest == 0] <- NA
expanse(full.forest, unit="ha", transform=TRUE, byValue=FALSE)

pre.fpct <- app(metric[[pre.idx]], fun=\(j) { mean(j >= forest.threshold) } ) 
  pre.forest <- mask(ifel(pre.fpct >= forest.pct, 1, 0), m)

post.fpct <- app(metric[[post.idx]], fun=\(j) { mean(j >= forest.threshold) } ) 
  post.forest <- ifel(post.fpct >= forest.pct, 1, 0)



fpct <- c(full.fpct, pre.fpct, post.fpct)
  names(fpct) <- c("full", "pre", "post")
forest <- c(full.forest, pre.forest, post.forest)
  names(forest) <- c("full", "pre", "post")
writeRaster(forest, file.path(results.dir, paste0("forest_mask_",mres,".tif")), 
  	        datatype = "INT2U", overwrite = TRUE)  
writeRaster(fpct, file.path(results.dir, paste0("forest_pct_",mres,".tif")), 
  	        datatype = "INT2U", overwrite = TRUE)  











x <- as.numeric(metric[154])	
s <- impute.loess(y=x, s = 0.2, smooth = TRUE)
dat <- data.frame(dates=dates, lai=s)	
	
ggplot(data=dat, aes(x=dates,y=lai, fill=lai))+
  geom_line()+
    geom_hline(yintercept=0)+
      geom_ribbon(data=dat,
              aes(ymin = lai, ymax = 0),
              fill="blue",alpha=0.5) +
	    theme_classic() + 
		  theme(rect = element_rect(fill = "transparent"))

theme(
  panel.background = element_rect(fill = "transparent", colour = NA_character_), 
  panel.grid.major = element_blank(), 
  panel.grid.minor = element_blank(), 
  plot.background = element_rect(fill = "transparent", colour = NA_character_), 
  legend.background = element_rect(fill = "transparent"),
  legend.box.background = element_rect(fill = "transparent"),
  legend.key = element_rect(fill = "transparent")
)

d <- density(dat$diff) 
d <- data.frame(x=d$x, x=d$y, Sign=sign(d$x))

ggplot(data=d, aes(x = x, y = dif, fill = dif))+
  geom_line()+
    #geom_hline(yintercept=0)+
      geom_ribbon(data = subset(d, Sign == -1),
              aes(xmin = dif, Xmax = 0),
              fill="red",alpha=0.5) +  
      geom_ribbon(data = subset(d, Sign == 1),
                  aes(xmin = 0, xmax = dif),
                  fill="blue",alpha=0.5)
				  
d <- density(dat$diff)
par(bg=NA)
shadeArea(d, shade.from = 0, shade.to = max(d$x), curve.col="grey",
          shade.col="darkseagreen")
shadeArea(d, shade.from = min(d$x), shade.to = -0.000001,
          curve.col="grey", shade.col="coral4", add = TRUE) 		  