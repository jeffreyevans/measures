suppressMessages(
  lapply(c("sf", "spatialEco", "terra","dplyr", "viridis", 
         "exactextractr", "ggplot2", "ggridges"), 
		 require, character.only = TRUE))

#*************************************************
# Set environment and model parameters 
metrics <- c("lai", "fcov")[2]

root = "C:/evans/Colombia/waterfunds"
setwd(root)
  dat.dir = file.path(getwd(), "data")
  mdl.dir <- file.path(getwd(), "model")

#**************************************
# read boundary and format interventions
vect.dat <- file.path(dat.dir, "CuencaVerde.gpkg") 
bdy <- st_read(vect.dat, "watersheds") |>
    st_union() |>
        st_sf() |> 
          st_cast(to="POLYGON")
    st_geometry(bdy) <- "geometry"

pa <- st_cast(st_read(vect.dat, layer="protected_areas"), 
              "POLYGON")
  st_geometry(pa) <- "geometry" 
    pa$intervention <- "Protected Area"
	  names(pa)[2] <- "name"
pa <- pa[,c("name", "intervention")]
st_geometry(pa) <- "geometry"

treatment <- st_read(vect.dat, "treatments") 
  treatment$intervention <- "Treatment"
  	names(treatment)[2] <- "name"
treatment <- treatment[,c("name", "intervention")]
st_geometry(treatment) <- "geometry"

controls <- st_read(vect.dat, "treatment_controls") 
  controls$intervention <- "Control"
    controls <- controls[,c("name", "intervention")]
st_geometry(controls) <- "geometry"

watersheds <- st_read(vect.dat, "watersheds") 

interventions <- unique(watersheds$label)
 
#**************************************
# read and classify effect sizes
es.raster <- merge(rast(file.path(mdl.dir, "causal", "intervention_fcov.tif")),
                   rast(file.path(mdl.dir, "causal", "Protected Area_fcov.tif")) )
  positive.es <- ifel(es.raster[[1]] > 0, 1, NA)
  negative.es <- ifel(es.raster[[1]] < 0, -1, NA)
  netural.es <- es.raster[[1]] 
    netural.es[!is.na(netural.es)] <- 0
  
  es.class <- merge(sprc(list(positive.es,negative.es,netural.es))) 
   levels(es.class) <- data.frame(Value=c(-1,0,1), 
         class = c("Negative effect", "No effect", "Positive effect"))

cols = c("darkviolet","bisque3","darkgreen")
plot(es.class, col=cols[c(1,3)])
  plot(st_geometry(watersheds),border="blue", add=TRUE)
  plot(st_geometry(pa),add=TRUE)

freq(es.class)
s <- cellSize(es.class, mask=TRUE, unit="km")
  es.class <- c(es.class, s)
esdf <- as.data.frame(es.class)
tapply(esdf$area, esdf$class, sum)

# order: Río Arma, Río Negro, Río Aburrá, Río Grande   
esc <- na.omit(extract(es.class, vect(watersheds)))
with(esc, tapply(area, list(ID, class), sum, na.rm=TRUE))

# PDF of effect size
d <- density(es[[1]])
plot(d[[1]], type="n", main="fractional cover effect size", xlab="PDF")
  polygon(d[[1]], col="blue")
    abline(v=0, col="grey")

#**************************************
# extract effect size by intervention(s)
# order: Río Arma, Río Negro, Río Aburrá, Río Grande   
units <- unique(watersheds$label)
v <- extract(es.raster[[1]], vect(watersheds))
  v <- lapply(unique(v$ID), \(i) {
    data.frame(intervention = unique(watersheds[i,]$label), es=as.numeric(v[v$ID == i,][,2]))
  })
  v <- do.call(rbind, v)
    v <- na.omit(v)
      v$intervention <- factor(v$intervention, levels=units)  

iord <- c("Río Grande", "Río Aburrá", "Río Negro", "Río Arma")
  v <- v[order(match(v$intervention, iord)),]
    v <- na.omit(v)

#**************************************
# Stacked quantile distributions
# using relative scale indicating observed effect sizes
qs <- seq(0.001, 0.95, 0.05)
# qs <- c(0.05, 0.10, 0.25, 0.5, 0.75, 0.90, 0.95)
ggplot(v, aes(x=es, y=intervention, fill=factor(..quantile..))) +
  stat_density_ridges(quantiles = qs, show.legend = FALSE,
                      quantile_lines = TRUE,
					  geom = "density_ridges_gradient") +
	scale_fill_viridis(discrete = TRUE, name = "Quantile", option = "plasma") +
	  #scale_fill_discrete(labels=c("-0.05-0", "0.10-0.25", "0.25-0.50", 
	  #                    "0.50-0.75", "0.75-0.90", "0.90-1")) +
		labs(fill="Quantiles") +
		  geom_vline(xintercept = 0, color = "black", size=1)
	

plot(tau.fcov[[2]], col=rev(grDevices::topo.colors(50)),
     alpha=0.6,legend=FALSE)
  plot(st_geometry(watersheds),add=TRUE)



#**************************************
# extract effect size for intervention(s)
i = interventions[1]
pa.int <- pa[pa$cat_manejo == i,]
  units <- unique(pa.int$nombre_asp)

v <- lapply(units, \(i) {
  na.omit(extract(es.raster[[1]], vect(pa.int[pa.int$nombre_asp == i,])))[,2]
})
names(v) <- units
  cts.idx <- which(unlist(lapply(v, \(i) length(i[!is.na(i)]))) > 0) 
    v <- v[cts.idx]
    for(i in 1:length(v)) {
      v[[i]] <- data.frame(ID=names(v)[i], es=v[[i]])
    }
v <- do.call(rbind, v)
  v$ID <- factor(v$ID)

#**************************************
# Stacked quantile distributions
# using fixed scale from, to of expected effect sizes (-1 - 1)

np <- floor(length(unique(v$ID)) / 4)
v.sub <- v[which(v$ID %in% unique(v$ID)[19:24]),]

ggplot(v, aes(x=es, y=intervention, fill=factor(..quantile..))) +
  stat_density_ridges(quantiles = c(0.1, 0.25, 0.5, 0.75, 0.9),
                      quantile_lines = TRUE, from=-1, to=1, 
					  geom = "density_ridges_gradient") +
	scale_fill_viridis(discrete = TRUE, name = "Quantile", option = "plasma") +
	  scale_fill_discrete(labels=c("0-0.10", "0.10-0.25", "0.25-0.50", 
	                      "0.50-0.75", "0.75-0.90", "0.90-1")) +
		labs(fill="Quantiles") 
	
	