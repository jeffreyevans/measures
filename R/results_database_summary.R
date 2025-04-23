# Effect size and trend results database query
#
# cell sizes in area 
#   250m = 0.0625 km, 6.25 ha
#   300m = 0.0899 km, 9 ha
#   500m = 0.2500 km  25 ha 
#
suppressMessages(
  lapply(c("sf", "spatialEco", "terra", "dplyr", "duckdb", "DBI",  
         "ggplot2", "rts", "fst", "parallel", "snow", "doSNOW", "xts"), 
		 require, character.only = TRUE))

#***************************************
#***************************************
# set up environment
mres = c("250m", "300m", "500m")[1]
m = c("lai", "fcov")[1]
if(mres == "250m") { 
    out.res = 250
    m = "lai"
	ha.scalar = 6.25 
} else if(mres == "500m") { 
    out.res = 500
    m = "fcov"
	ha.scalar = 25	
} else if(mres == "300m") {  
    out.res = 300 
	ha.scalar = 9
 }

idx = 2  # indicates country 
cty <- c("bhutan", "costa_rica", "canada", "colombia", "peru", "brazil")[idx]

root = file.path("C:/evans/PFP", cty)
  mdl.dir = file.path(root, "model", paste0("model",mres))
  dat.dir = file.path(root, "data")
  setwd(mdl.dir)

#******************************************************
# Read polygon data and build PA/PFP predicates
bdy <- st_read(file.path(dat.dir, paste0(cty, ".gpkg")), "boundary")
if(any(cty %in% c("colombia", "peru", "brazil"))){ 	  
  pa <- st_read(file.path(dat.dir, paste0(cty, ".gpkg")), "protected_areas")
    st_geometry(pa) <- "geometry" 
      if(cty == "brazil") pa <- pa[-which(pa$DESIG %in% "Indigenous Area"),]
  	pa$type <- "protected Area"
  pfp <- st_read(file.path(dat.dir, paste0(cty, ".gpkg")), "PFP")
    st_geometry(pfp) <- "geometry"
  	pfp$IUCN <- "Not Applicable"
	pfp$type <- "PFP"
  pfo.wo.pa <- st_read(file.path(dat.dir, paste0(cty, ".gpkg")), "pfp_wo_pa") 
  pa.with.pfp <- pa[which(lengths(st_intersects(pa, pfp)) > 0),]  
} else {
  pa <- st_read(file.path(dat.dir, paste0(cty, ".gpkg")), "protected_areas")
    st_geometry(pa) <- "geometry"
}

# create IUCN aggregrated polygons (collapse Ia & Ib)
pa$IUCN[which(pa$IUCN %in% "Ia" | pa$IUCN %in% "Ib")] <- "I"
  iucn <- pa %>%
    dplyr::group_by(IUCN) %>% 
      summarize(geometry = st_union(geometry))
plot(iucn["IUCN"], border=NA)

#******************************************************
# read lai and fractional cover
# Date ranges
date.break = c("2017-01-01", "2010-01-01", "2006-01-01", "2020-01-01", "2019-01-01", "2014-01-01")[idx]
ts.date <- list(c(2017, 1, 1), c(2010, 1, 1), c(2006, 1, 1), c(2020, 1, 1), c(2019, 1, 1), c(2014, 1, 1))
  names(ts.date) <- c("bhutan", "costa_rica", "canada", "colombia", "peru", "brazil") 

# read lai timeseries and set pre/post intervention dates
metric <- rast(file.path(dat.dir, "LAI_250m_trend.tif"))
  dates <- as.Date(names(metric))
    pre.idx <- which(dates < date.break)
    post.idx <- which(dates >= date.break)

#******************************************************
# Amount of forest in study period
msk <- rast(file.path(dat.dir, paste0("mask", mres, ".tif")))
f <- rast(file.path(dat.dir, paste0("forest_", mres, ".tif")))
  f[f == 0] <- NA
    a <- mask(cellSize(f, unit="ha"), f)
global(a, fun="sum", na.rm = TRUE)   # km of forest
global(a, fun="mean", na.rm = TRUE)  # mean cell size in km

# Forest change
baseline <- mean(range(metric[[grep("2000", dates)]], na.rm=TRUE))
current <- mean(range(metric[[grep("2021", dates)]], na.rm=TRUE))

# Country wide
change <- as.data.frame(mask(current - baseline, f))
  names(change) <- "change"
    change$change_class <- ifelse(change$change >= -0.5 & change$change <= 0.5, 0,
                              ifelse(change$change > -0.5, 1,
                                ifelse(change$change < -0.5, -1, NA)))
    chg.cts <- as.data.frame(table(change$change_class))
      names(chg.cts) <- c("change", "counts")
        chg.cts$pct <- chg.cts$counts / sum(chg.cts$counts)
    	chg.cts$ha <- chg.cts$counts * ha.scalar

# Protected areas
change <- as.data.frame(mask(current - baseline, pa))
  names(change) <- "change"
    change$change_class <- ifelse(change$change >= -0.5 & change$change <= 0.5, 0,
                             ifelse(change$change > -0.5, 1,
                               ifelse(change$change < -0.5, -1, NA)))
    chg.cts <- as.data.frame(table(change$change_class))
      names(chg.cts) <- c("change", "counts")
        chg.cts$pct <- chg.cts$counts / sum(chg.cts$counts)
    	chg.cts$ha <- chg.cts$counts * ha.scalar

length(which(change$change > 2)) * ha.scalar
length(which(change$change < -2)) * ha.scalar
	
	
#******************************************************
# Read results database
results.db <- file.path(mdl.dir, paste0(m, "_", mres, "_results.duckdb"))   
  con <- dbConnect(duckdb::duckdb(dbdir = results.db), read_only = FALSE)
    db.table <- dbListTables(con)
    dat <- dbGetQuery(con, paste0("SELECT * FROM ", db.table)) 
  dbDisconnect(con, shutdown = TRUE)

control <- dat[dat$NAME == "control",]
  control <- control[,-c(6:8,15:57)]
treatment <- dat[dat$NAME != "control",]
  treatment <- treatment[,-c(6:8,15:57)]
remove(dat) 

# for IUCN or elu analysis, subset here
# will have to reread above data for different
# analysis

# treatment <- treatment[grep(paste(c("I", "Ia", "Ib"), collapse="|"), treatment$IUCN),]
# treatment <- treatment[grep(paste(c("I", "Ia", "Ib", "II"), collapse="|"), treatment$IUCN),]

#################################
####  treatment query
es.netural <- list(current=which(treatment$effect_size_current > -0.01 & treatment$effect_size_current < 0.01), 
                   current.pct=which(treatment$effect_size_cpct > -0.01 & treatment$effect_size_cpct < 0.01),
                   tsa=which(treatment$effect_size_tsa > -0.01 & treatment$effect_size_tsa < 0.01))

es.pos <- list(current=which(treatment$effect_size_current >= 0.01), 
               current.pct=which(treatment$effect_size_cpct >= 0.01),
               tsa=which(treatment$effect_size_tsa >= 0.01))
			   
es.neg <- list(current=which(treatment$effect_size_current <= -0.01), 
               current.pct=which(treatment$effect_size_cpct <= -0.01),
               tsa=which(treatment$effect_size_tsa <= -0.01)) 
			   
es.high.neg <- list(current=which(treatment$effect_size_current <= -0.5), 
                    current.pct=which(treatment$effect_size_cpct <= -0.5),
                     tsa=which(treatment$effect_size_tsa <= -0.5))	
					 
es.high.pos <- list(current=which(treatment$effect_size_current >= 0.5), 
                    current.pct=which(treatment$effect_size_cpct >= 0.5),
                     tsa=which(treatment$effect_size_tsa >= 0.5))
					 
es.cts <- data.frame(class=c("neg", "neg_high", "netural", "pos", "pos_high"),
                     rbind(unlist(lapply(es.neg, length)),
                     unlist(lapply(es.high.neg, length)),
					 unlist(lapply(es.netural, length)), 
                     unlist(lapply(es.pos, length)),
                     unlist(lapply(es.high.pos, length))))

n = length(which(!is.na(treatment$effect_size_current)))
es.pct <- data.frame(class=es.cts$class, es.cts[,-1] / n) 
es.ha <- data.frame(class=es.cts$class, es.cts[,-1] * ha.scalar)

#****************************************************
#****************************************************
# Trend counts within positive treatmets
trend.names <- c("pre_high_neg_trend", "pre_neg_trend", 
                 "pre_netural_trend", "pre_pos_trend", 
                 "pre_pos_high_trend", "post_high_neg_trend", 
                 "post_neg_trend", "post_netural_trend", 
                 "post_pos_trend", "post_pos_high_trend",
                 "neg_trend_high_effect", "high_pos_trend_high_effect")

#### Current tau
tau.pre <- treatment$tau_pre[es.pos[["current"]]]
tau.post <- treatment$tau_post[es.pos[["current"]]]
curr.trend <- na.omit(data.frame(intervention=c(rep("pre",length(tau.pre)), 
                         rep("post",length(tau.post))),
                         tau=c(tau.pre,tau.post)))
post.tau.es50 <- treatment$tau_post[which(treatment$effect_size_current >= 0.50)]
n = length(tau.post) 
  curr.trend.cts <- 
    c(length(which(curr.trend[curr.trend$intervention == "pre",]$tau <= -0.30)),
	
      length(which(curr.trend[curr.trend$intervention == "pre",]$tau <= -0.05)),
	  
      length(which(curr.trend[curr.trend$intervention == "pre",]$tau > -0.05 & 
                   curr.trend[curr.trend$intervention == "pre",]$tau < 0.05)), 
				   
      length(which(curr.trend[curr.trend$intervention == "pre",]$tau >= 0.05)),
	  
      length(which(curr.trend[curr.trend$intervention == "pre",]$tau >= 0.30)),
	  
	  
      length(which(curr.trend[curr.trend$intervention == "post",]$tau <= -0.30)),
      length(which(curr.trend[curr.trend$intervention == "post",]$tau <= -0.05)),
      length(which(curr.trend[curr.trend$intervention == "post",]$tau > -0.05 & 
                   curr.trend[curr.trend$intervention == "post",]$tau < 0.05)), 
      length(which(curr.trend[curr.trend$intervention == "post",]$tau >= 0.05)),
      length(which(curr.trend[curr.trend$intervention == "post",]$tau >= 0.50)),
      length(which(post.tau.es50 < -0.01)),
      length(which(post.tau.es50 > 0.3)) )
names(curr.trend.cts) <- trend.names 

#### Current-percent tau
tau.pre <- treatment$tau_pre[es.pos[["current.pct"]]]
tau.post <- treatment$tau_post[es.pos[["current.pct"]]]
cpct.trend <- na.omit(data.frame(intervention=c(rep("pre",length(tau.pre)), 
                         rep("post",length(tau.post))),
                         tau=c(tau.pre,tau.post)))
post.tau.es50 <- treatment$tau_post[which(treatment$effect_size_cpct >= 0.50)]
n = length(tau.post) 
  cpct.trend.cts <-
    c(length(which(cpct.trend[cpct.trend$intervention == "pre",]$tau <= -0.30)),
      length(which(cpct.trend[cpct.trend$intervention == "pre",]$tau <= -0.05)),
      length(which(cpct.trend[cpct.trend$intervention == "pre",]$tau > -0.05 & 
                   cpct.trend[cpct.trend$intervention == "pre",]$tau < 0.05)), 
      length(which(cpct.trend[cpct.trend$intervention == "pre",]$tau >= 0.05)),
      length(which(cpct.trend[cpct.trend$intervention == "pre",]$tau >= 0.30)),
      length(which(cpct.trend[cpct.trend$intervention == "post",]$tau <= -0.30)),
      length(which(cpct.trend[cpct.trend$intervention == "post",]$tau <= -0.05)),
      length(which(cpct.trend[cpct.trend$intervention == "post",]$tau > -0.05 & 
                   cpct.trend[cpct.trend$intervention == "post",]$tau < 0.05)), 
      length(which(cpct.trend[cpct.trend$intervention == "post",]$tau >= 0.05)),
      length(which(cpct.trend[cpct.trend$intervention == "post",]$tau >= 0.50)),
      length(which(post.tau.es50 < -0.01)),                             
      length(which(post.tau.es50 > 0.3)) )
names(cpct.trend.cts) <- trend.names 

#### temporal-volume tau
tau.pre <- treatment$tau_pre[es.pos[["tsa"]]]
tau.post <- treatment$tau_post[es.pos[["tsa"]]]
tsa.trend <- na.omit(data.frame(intervention=c(rep("pre",length(tau.pre)), 
                         rep("post",length(tau.post))),
                         tau=c(tau.pre,tau.post)))
post.tau.es50 <- treatment$tau_post[which(treatment$effect_size_tsa >= 0.50)]
n = length(tau.post) 						 
  tsa.trend.cts <-
    c(length(which(tsa.trend[tsa.trend$intervention == "pre",]$tau  <= -0.30)),
      length(which(tsa.trend[tsa.trend$intervention == "pre",]$tau  <= -0.05)),
      length(which(tsa.trend[tsa.trend$intervention == "pre",]$tau  > -0.05 & 
                   tsa.trend[tsa.trend$intervention == "pre",]$tau  < 0.05)), 
      length(which(tsa.trend[tsa.trend$intervention == "pre",]$tau  >= 0.05)),
      length(which(tsa.trend[tsa.trend$intervention == "pre",]$tau  >= 0.30)),
      length(which(tsa.trend[tsa.trend$intervention == "post",]$tau <= -0.30)),
      length(which(tsa.trend[tsa.trend$intervention == "post",]$tau <= -0.05)),
      length(which(tsa.trend[tsa.trend$intervention == "post",]$tau > -0.05 & 
                   tsa.trend[tsa.trend$intervention == "post",]$tau < 0.05)), 
      length(which(tsa.trend[tsa.trend$intervention == "post",]$tau >= 0.05)),
      length(which(tsa.trend[tsa.trend$intervention == "post",]$tau >= 0.50)),
      length(which(post.tau.es50 < -0.01)),
      length(which(post.tau.es50 > 0.3)) )  
names(tsa.trend.cts) <- trend.names   


   pre.cts <- setNames(data.frame(cbind(curr.trend.cts[grep("pre", names(curr.trend.cts))], 
	  cpct.trend.cts[grep("pre", names(cpct.trend.cts))], 
	  tsa.trend.cts[grep("pre", names(tsa.trend.cts))])),
	  c("current", "current percent", "biomass volume"))
        pre.cts <- data.frame(period="pre intervention", trend = rownames(pre.cts), pre.cts) 
   post.cts <- setNames(data.frame(cbind(curr.trend.cts[grep("post", names(curr.trend.cts))], 
	  cpct.trend.cts[grep("pre", names(cpct.trend.cts))], 
	  tsa.trend.cts[grep("pre", names(tsa.trend.cts))])),
	  c("current", "current percent", "biomass volume")) 
    post.cts <- data.frame(period="post intervention", trend = rownames(post.cts), post.cts) 

trend.cts <- rbind(pre.cts, post.cts)
  rownames(trend.cts) <- 1:nrow(trend.cts)

curr.den <- ggplot(curr.trend, aes(x=tau, fill=intervention)) +
  geom_density(alpha=0.3) +
     labs(title="Trends in negative current", x="trend", y = "density")
cpct.den <- ggplot(cpct.trend, aes(x=tau, fill=intervention)) +
  geom_density(alpha=0.3) +
     labs(title="Trends in negative current-pct", x="trend", y = "density")
tsa.den <- ggplot(tsa.trend, aes(x=tau, fill=intervention)) +
  geom_density(alpha=0.3) +
     labs(title="Trends in negative biomass volume", x="trend", y = "density")  

ggpubr::ggarrange(curr.den, cpct.den, tsa.den)

#################################
#### negative treatmets
# Current tau
tau.pre <- treatment$tau_pre[es.neg[["current"]]]
tau.post <- treatment$tau_post[es.neg[["current"]]]
curr.trend <- data.frame(intervention=c(rep("pre",length(tau.pre)), 
                         rep("post",length(tau.post))),
                         tau=c(tau.pre,tau.post))
  length(which(curr.trend$pre <= -0.05)) * ha.scalar
  length(which(curr.trend$pre <= -0.05)) / n
  
  length(which(curr.trend$pre > -0.05 & curr.trend$pre < 0.05)) * ha.scalar
  length(which(curr.trend$pre > -0.05 & curr.trend$pre < 0.05)) / n  
  
  length(which(curr.trend$pre >= 0.05)) * ha.scalar
  length(which(curr.trend$pre >= 0.05)) / n

# Current-percent tau
tau.pre <- treatment$tau_pre[es.neg[["current.pct"]]]
tau.post <- treatment$tau_post[es.neg[["current.pct"]]]
cpct.trend <- na.omit(data.frame(pre=tau.pre, post=tau.post))
  length(which(cpct.trend$pre <= -0.05)) * ha.scalar
  length(which(cpct.trend$pre > -0.05 & cpct.trend$pre < 0.05)) * ha.scalar
  length(which(cpct.trend$pre >= 0.05)) * ha.scalar

# temporal-volume tau
tau.pre <- treatment$tau_pre[es.neg[["tsa"]]]
tau.post <- treatment$tau_post[es.neg[["tsa"]]]
tsa.trend <- na.omit(data.frame(pre=tau.pre, post=tau.post))
  length(which(tsa.trend$pre <= -0.05)) * ha.scalar
  length(which(tsa.trend$pre > -0.05 & tsa.trend$pre < 0.05)) * ha.scalar
  length(which(tsa.trend$pre >= 0.05)) * ha.scalar
