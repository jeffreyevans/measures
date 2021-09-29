#  eu - sp polygons of experimental units (interventions)
#      IND002.3 - Community Forestry Areas  
#      IND002.4 - Community Initiative Berau Barat (KPHP)
#      IND002.5 - Community Initiative Wehea-Kelay Ecosystem Essential Area
#      Protected Areas - National parks   
################################################################

load("C:/evans/measures/Indonesia/CIA_PixelLevel_model.RData")
   
library(sp)
library(raster)
library(spatialEco)
library(rgdal)
library(sf)
library(lwgeom)
library(RColorBrewer)

setwd("C:/evans/measures/Indonesia")
  fdir = file.path(getwd(), "figures/plots")

#### boundary data used in analysis 
bdy <- readOGR(file.path(getwd(), "boundaries"),"KalimantanTimur")

# Fix invalid topography in experimental unit polygons   
# and explode features 
st_is_valid(as(eu, "sf"), reason = TRUE)
  eu.sf <- as(eu, "sf") %>% 
    lwgeom::st_make_valid() %>% 
	  sf::st_cast("POLYGON")
eu <- as(eu.sf, "Spatial")

yr.idx <- grep("2008", lai_dates)
ldf <- point.in.poly(lai.all, eu)
  ldf <- sp.na.omit(ldf)
kendall.tau <- as.data.frame(t(apply(ldf@data[1:576][yr.idx], 
                             MARGIN=1, FUN=kendall)))
	names(kendall.tau)[4] <- "p"
  ldf@data <- data.frame(ldf@data[,577:580], kendall.tau) 


trend <- as(ldf,"sf")

es.pal <- rev(brewer.pal(10, "Spectral"))
  plot(trend["tau"], pal = es.pal, pch=20, cex=0.5, 
       breaks = seq(-1,1,0.2), axes = TRUE, 
	   border=NA,, key.pos = 4)

###################################
# TAU
### All
n=dim(ldf)[1]
all.up = dim(ldf[ldf$tau >= 0.20,])[1]
  all.up = all.up / n
all.down = dim(ldf[ldf$tau <= -0.2,])[1]
  all.down = all.down / n
all.stable = 1 - (all.up+all.down) 

### IND002.3 - Community Forestry Areas 
sf.trend <- ldf@data[ldf@data$PROJECT_ID == "IND002.3",]
n=dim(sf.trend)[1]
  sf.up = dim(sf.trend[sf.trend$tau >= 0.20,])[1]
    sf.up = sf.up / n
  sf.down = dim(sf.trend[sf.trend$tau <= -0.2,])[1]
    sf.down = sf.down / n
  sf.stable = sf.stable = 1 - (sf.up+sf.down) 

### IND002.4 - Community Initiative Berau Barat (KPHP)
ma.trend <- ldf@data[ldf@data$PROJECT_ID == "IND002.4",]
n=dim(ma.trend)[1]
  ma.up = dim(ma.trend[ma.trend$tau >= 0.20,])[1]
    ma.up = ma.up / n
  ma.down = dim(ma.trend[ma.trend$tau <= -0.2,])[1]
    ma.down = ma.down / n
  ma.stable = ma.stable = 1 - (ma.up+ma.down) 

### IND002.5 - Community Initiative Wehea-Kelay Ecosystem Essential Area
ea.trend <- ldf@data[ldf@data$PROJECT_ID == "IND002.5",]
n=dim(ea.trend)[1]
  ea.up = dim(ea.trend[ea.trend$tau >= 0.20,])[1]
    ea.up = ea.up / n
  ea.down = dim(ea.trend[ea.trend$tau <= -0.2,])[1]
    ea.down = ea.down / n
  ea.stable = ea.stable = 1 - (ea.up+ea.down) 
  
### Protected Areas - National parks
pa.trend <- ldf@data[ldf@data$PROJECT_ID == "Protected Areas",]
n=dim(pa.trend)[1]
  pa.up = dim(pa.trend[pa.trend$tau >= 0.20,])[1]
    pa.up = pa.up / n
  pa.down = dim(pa.trend[pa.trend$tau <= -0.20,])[1]
    pa.down = pa.down / n
  pa.stable = pa.stable = 1 - (pa.up+pa.down) 

tau.pct <- rbind(all.up = all.up, 
                all.stable = all.stable,
                all.down = all.down,
				sf.up = sf.up, 
				sf.stable = sf.stable,
                sf.down = sf.down,
                ma.up = ma.up, 
				ma.stable = ma.stable,
                ma.down = ma.down,	
                ea.up = ea.up, 
				ea.stable = ea.stable,
                ea.down = ea.down,
                pa.up = pa.up, 
				pa.stable = pa.stable,
                pa.down = pa.down)
	 
###################################
# slope

### All
n=dim(ldf)[1]
all.up = dim(ldf[ldf$slope > 0,])[1]
  all.up = all.up / n
all.down = dim(ldf[ldf$slope < 0,])[1]
  all.down = all.down / n
all.stable = 1 - (all.up+all.down) 

### IND002.3 - Community Forestry Areas 
sf.trend <- ldf@data[ldf@data$PROJECT_ID == "IND002.3",]
n=dim(sf.trend)[1]
  sf.up = dim(sf.trend[sf.trend$slope > 0,])[1]
    sf.up = sf.up / n
  sf.down = dim(sf.trend[sf.trend$slope < 0,])[1]
    sf.down = sf.down / n
  sf.stable = sf.stable = 1 - (sf.up+sf.down) 

### IND002.4 - Community Initiative Berau Barat (KPHP)
ma.trend <- ldf@data[ldf@data$PROJECT_ID == "IND002.4",]
n=dim(ma.trend)[1]
  ma.up = dim(ma.trend[ma.trend$slope > 0,])[1]
    ma.up = ma.up / n
  ma.down = dim(ma.trend[ma.trend$slope < 0,])[1]
    ma.down = ma.down / n
  ma.stable = ma.stable = 1 - (ma.up+ma.down) 

### IND002.5 - Community Initiative Wehea-Kelay Ecosystem Essential Area
ea.trend <- ldf@data[ldf@data$PROJECT_ID == "IND002.5",]
n=dim(ea.trend)[1]
  ea.up = dim(ea.trend[ea.trend$slope > 0,])[1]
    ea.up = ea.up / n
  ea.down = dim(ea.trend[ea.trend$slope < 0,])[1]
    ea.down = ea.down / n
  ea.stable = ea.stable = 1 - (ea.up+ea.down) 
  
### Protected Areas - National parks
pa.trend <- ldf@data[ldf@data$PROJECT_ID == "Protected Areas",]
n=dim(pa.trend)[1]
  pa.up = dim(pa.trend[pa.trend$slope > 0,])[1]
    pa.up = pa.up / n
  pa.down = dim(pa.trend[pa.trend$slope < 0,])[1]
    pa.down = pa.down / n
  pa.stable = pa.stable = 1 - (pa.up+pa.down) 

slp.pct <- rbind(all.up = all.up, all.stable = all.stable,
                      all.down = all.down,
					  sf.up = sf.up, sf.stable = sf.stable,
                      sf.down = all.down,
                      ma.up = ma.up, ma.stable = ma.stable,
                      ma.down = all.down,					
                      ea.up = ea.up, ea.stable = ea.stable,
                      ea.down = all.down,
                      pa.up = pa.up, pa.stable = pa.stable,
                      pa.down = all.down)
	                             
trend.pct <- data.frame(type=rownames(tau.pct), 
                        tau=tau.pct[,1], 
			            slope=slp.pct[,1])

