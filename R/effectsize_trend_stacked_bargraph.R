library(ggplot2)

setwd("C:/evans/PFP/results/tables")

# effect size data
dat <- read.csv("lai_interventions_effectsize_barplots.csv")
  dat.chg <- dat[which(dat$metric %in% "change"),][,-2]
  dat.curr <- dat[which(dat$metric %in% "current"),][,-2]
  dat.vol <- dat[which(dat$metric %in% "tsa"),][,-2]

# trend data
dat.trend <- read.csv("lai_interventions_trend_barplots.csv")
  dat.trend <- dat.trend[dat.trend$trend != "small",]
  dat.trend$trend <- factor(dat.trend$trend)
    dat.pre <- dat.trend[,-4]
      dat.pre$period <- "pre"
      names(dat.pre)[3] <- "percent"
    dat.post <- dat.trend[,-3]
      dat.post$period <- "post"
      names(dat.post)[3] <- "percent"
    dat.trend <- rbind(dat.pre, dat.post)
      dat.trend$period <- factor(dat.trend$period, levels = c("pre","post"))

pdf("effect_size_trend_bargraphs.pdf", height=8.5, width=11)
  chg = ggplot(dat.chg, aes(x = country, y = percent, fill = effect_size)) +  
        geom_bar(data = dat.chg[dat.chg$effect_size != "small",], 
                 position = "stack", stat = "identity") + 
        scale_fill_manual(values=c("darkred", "darkseagreen4")) +
        geom_col(data = dat.chg[dat.chg$effect_size == "small", ], 
		         aes(group = effect_size), fill = alpha("grey", 0.45), 
				 width = 0.5, position = position_dodge(width = 0.9)) +
        theme_minimal() +
    	ggtitle("Absolute Change (positive, negative and small values around zero)") +
    	xlab("") + ylab("percent") + guides(fill=guide_legend(title="")) 
  print(chg)

  curr = ggplot(dat.curr, aes(x = country, y = percent, fill = effect_size)) +  
        geom_bar(data = dat.curr[dat.curr$effect_size != "small",], 
                 position = "stack", stat = "identity") + 
        scale_fill_manual(values=c("darkred", "darkseagreen4")) +
        geom_col(data = dat.curr[dat.curr$effect_size == "small", ], 
		         aes(group = effect_size), fill = alpha("grey", 0.45), 
				 width = 0.5, position = position_dodge(width = 0.9)) +
        theme_minimal() +
    	ggtitle("currrent effect size (positive, negative and small values around zero)") +
    	xlab("") + ylab("percent") + guides(fill=guide_legend(title=""))  
  print(curr)		

  vol = ggplot(dat.vol, aes(x = country, y = percent, fill = effect_size)) +  
        geom_bar(data = dat.vol[dat.vol$effect_size != "small",], 
                 position = "stack", stat = "identity") + 
        scale_fill_manual(values=c("darkred", "darkseagreen4")) +
        geom_col(data = dat.vol[dat.vol$effect_size == "small", ], 
		         aes(group = effect_size), fill = alpha("grey", 0.45), 
				 width = 0.5, position = position_dodge(width = 0.9)) +
        theme_minimal() +
    	ggtitle("Temporal volume effect size (positive, negative and small values around zero)") +
    	xlab("") + ylab("percent") + guides(fill=guide_legend(title=""))  
  print(vol)

  trend <- ggplot(dat.trend, aes(x = period, y = percent, fill = trend)) +
    geom_bar(position = "stack", stat = "identity") +
    scale_fill_manual(values=c("darkred", "darkseagreen4")) +		      
    facet_grid(~ country, switch = "x") +
    theme(strip.placement = "outside",
          strip.background = element_rect(fill = NA, color = "white"),
          panel.spacing = unit(-0.01,"cm")) +
  		theme_minimal() +
  		ggtitle("pre and post intervention percent positive/negative LAI trends (interventions)") +
  		xlab("") + ylab("percent")
  print(trend)  
dev.off()
