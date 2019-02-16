### Created by SHC May 2018
source("analyze_data/egg_analysis_build_dataframe.R")
source("analyze_data/egg_analysis_read_trees.R")

### This code plots a comparison between egg features based on oviposition ecology
analysis_name <- "plotting_ecology"
class_flag <- "relaxed"

ecology_arg <- "internal"
source("analyze_data/egg_analysis_parasitoid.R")
plotting_data1 <- egg_eco_data %>% filter(!(group %in% c("Apterygota","Psocodea")))
ecology_arg <- "aquatic"
source("analyze_data/egg_analysis_aquatic.R")

### Set up data and plot parameters
plotting_data2 <- egg_eco_data %>% filter(!(group %in% c("Apterygota","Psocodea")))
vol_breaks <- c(10^-6,10^-5,10^-4,10^-3,10^-2,10^-1,1,10,100,1000,10000,100000,1000000)
ar_breaks <- c(0.25,0.5,1,2,4,8,16)
asym_breaks <- c(0,0.04,.16,.36,.64,1)
asym_breaks_untransformed <- c(0,.25,.5,.75,1)
curv_breaks <- deg2rad(c(0.0,7.2,28.8,64.8,115.2,180.0))
curv_breaks_untransformed <- deg2rad(c(0,60,120,180))

eco_vol_strip_plot <- ggplot(plotting_data2,aes(x = factor(group,levels=rev(group_levels)), y = vol)) + 
	geom_jitter(data = plotting_data2,pch=NA) + 
	geom_jitter(data = plotting_data2 %>% filter(eco_regime == "ancestral"),pch=16,color = "#d4d4d4",) + 
	geom_jitter(data = plotting_data2 %>% filter(eco_regime == "derived"),pch=16,color = "#2056CE") + 
	geom_jitter(data = plotting_data1 %>% filter(eco_regime == "derived"),pch=16,color = "#c55c15") + 
	scale_y_log10(breaks = vol_breaks) + 
	coord_flip() + 
	theme(legend.position = "none")


pdf(file=paste("eco_strip_vol_plot.pdf",sep="_"),width=9,height=6)
	print(eco_vol_strip_plot)
dev.off()

eco_ar_strip_plot <- ggplot(plotting_data2,aes(x = factor(group,levels=rev(group_levels)), y = ar)) + 
	geom_jitter(data = plotting_data2,pch=NA) + 
	geom_jitter(data = plotting_data2 %>% filter(eco_regime == "ancestral"),pch=16,color = "#d4d4d4",) + 
	geom_jitter(data = plotting_data2 %>% filter(eco_regime == "derived"),pch=16,color = "#2056CE") + 
	geom_jitter(data = plotting_data1 %>% filter(eco_regime == "derived"),pch=16,color = "#c55c15") + 
	scale_y_log10(breaks = ar_breaks) + 
	coord_flip() + 
	theme(legend.position = "none")

pdf(file=paste("eco_strip_ar_plot.pdf",sep="_"),width=9,height=6)
	print(eco_ar_strip_plot)
dev.off()

eco_curv_strip_plot <- ggplot(plotting_data2,aes(x = factor(group,levels=rev(group_levels)), y = curv)) + 
	geom_jitter(data = plotting_data2,pch=NA) + 
	geom_jitter(data = plotting_data2 %>% filter(eco_regime == "ancestral"),pch=16,color = "#d4d4d4",) + 
	geom_jitter(data = plotting_data2 %>% filter(eco_regime == "derived"),pch=16,color = "#2056CE") + 
	geom_jitter(data = plotting_data1 %>% filter(eco_regime == "derived"),pch=16,color = "#c55c15") + 
	scale_y_continuous(breaks = curv_breaks_untransformed) + 
	coord_flip() + 
	theme(legend.position = "none")

pdf(file=paste("eco_strip_curv_plot.pdf",sep="_"),width=9,height=6)
	print(eco_curv_strip_plot)
dev.off()

eco_asym_strip_plot <- ggplot(plotting_data2,aes(x = factor(group,levels=rev(group_levels)), y = asym)) + 
	geom_jitter(data = plotting_data2,pch=NA) + 
	geom_jitter(data = plotting_data2 %>% filter(eco_regime == "ancestral"),pch=16,color = "#d4d4d4",) + 
	geom_jitter(data = plotting_data2 %>% filter(eco_regime == "derived"),pch=16,color = "#2056CE") + 
	geom_jitter(data = plotting_data1 %>% filter(eco_regime == "derived"),pch=16,color = "#c55c15") + 
	scale_y_continuous(breaks = asym_breaks_untransformed) + 
	coord_flip() + 
	theme(legend.position = "none")

pdf(file=paste("eco_strip_asym_plot.pdf",sep="_"),width=9,height=6)
	print(eco_asym_strip_plot)
dev.off()