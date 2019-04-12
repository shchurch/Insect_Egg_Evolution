### Created by SHC April 2018

### This code plots the morphospcae of the insect egg database
### It also compares the insect morphospace to the published avian egg morphospace

source("analyze_data/egg_analysis_build_dataframe.R")
source("analyze_data/egg_analysis_convert_stoddard.R")

# Set morphospace tick marks
vol_breaks <- c(10^-6,10^-5,10^-4,10^-3,10^-2,10^-1,1,10,100,1000,10000,100000,1000000)
ar_breaks <- c(0.25,0.5,1,2,4,8,16)
asym_breaks <- c(0,0.04,.16,.36,.64,1)
asym_breaks_untransformed <- c(0,.25,.5,.75,1)
curv_breaks <- deg2rad(c(0.0,7.2,28.8,64.8,115.2,180.0))
curv_breaks_untransformed <- deg2rad(c(0,50,100,150,200))

# Plot polyembryony on morphospace
polyembryonic <- read.delim("analyze_data/polyembryonic_taxa.tsv",header=T,stringsAsFactors=F)

# set default to monoembryonic
egg_database$polyembryonic <- "monoembryonic"
# code in polyembryony
egg_database[which(egg_database$genus %in% polyembryonic$genus),]$polyembryonic <- "genus_has_polyembryony"
egg_database[which(egg_database$name %in% polyembryonic$name),]$polyembryonic <- "polyembryonic"

polyembryonic_plot <- ggplot(egg_database,aes(x = ar, y = vol,color = polyembryonic)) + 
	geom_point(color = "dark grey",alpha=0.9,pch=16) + 
	geom_point(data = egg_database %>% filter(!(polyembryonic == "monoembryonic"))) + 
	scale_x_log10(breaks = ar_breaks) + 
	scale_y_log10(breaks = vol_breaks,expand=c(0,0.2)) +
	scale_color_manual(values = c("#2056CE","#c55c15")) + 
	theme(legend.position = "none")

# Plot ootheca on morphospace
source("analyze_data/egg_analysis_ootheca.R")
egg_eco_data <- egg_eco_data %>% mutate(ootheca = ifelse(eco_regime == "ancestral","no ootheca","ootheca"))

ootheca_plot <- ggplot(egg_eco_data,aes(x = ar, y = vol,color = ootheca)) + 
	geom_point(color = "dark grey",alpha=0.9,pch=16) + 
	geom_point(data = egg_eco_data %>% filter(!(ootheca == "no ootheca"))) + 
	scale_x_log10(breaks = ar_breaks) + 
	scale_y_log10(breaks = vol_breaks,expand=c(0,0.2)) +
	scale_color_manual(values = c("#2056CE","#c55c15")) + 
	theme(legend.position = "none")

# Build pairwise morphospace plots with two parameters
morpho_volume_aspect_ratio <- ggplot(egg_database,aes(x = ar, y = vol,color = group)) + 
	geom_point(size=1.75,alpha = 0.65,pch=16) +	
	scale_x_log10(breaks = ar_breaks) + 
	scale_y_log10(breaks = vol_breaks,expand=c(0,0.2)) +
	theme(legend.position="none") + 
	scale_color_manual(values = mrk) 

# Pairwise plots by subgroup
submrk <-  c("#D454E2","#2056CE","#2E9B00","#22CED5","#BE0000","#E87425","#878787","#3A175F","#F2C406","#000000")

subgroup1 <- egg_database %>% 
	filter(group == "Apterygota") %>% 
	mutate(group = order) %>% 
	filter(group != "")

subgroup2 <- egg_database %>% 
	filter(group == "Palaeoptera") %>% 
	mutate(group = ifelse(order == "Ephemeroptera",
		"Ephemeroptera",suborder)) %>% 
	filter(group != "")

subgroup3 <- egg_database %>% 
	filter(group == "Polyneoptera") %>% 
	mutate(group = order) %>% 
	filter(group != "")

subgroup4 <- egg_database %>% 
	filter(group == "Condylognatha") %>% 
	mutate(group = suborder) %>% mutate(group = 
	ifelse(superfamily %in% c("Cercopoidea","Fulgoroidea","Membracoidea"),
	"Auchenorryncha",ifelse(order == "Thysanoptera",
		order,suborder))) %>% 
	filter(group != "") 

subgroup5 <- egg_database %>% 
	filter(group == "Hymenoptera") %>% 
	mutate(group = ifelse(
		superfamily %in% c("Tenthredinoidea", "Siricoidea", "Cephoidea", "Xiphidrioidea", 
		"Xyeloidea", "Orussoidea","Pamphilioidea","Xiphydrioidea"),
			"Symphyta",ifelse(
			superfamily %in% c("Ceraphronoidea","Chalcidoidea","Cynipoidea","Evanioidea",
			"Ichneumonoidea","Platygastroidea","Proctotrupoidea","Trigonaloidea"),
				"Parasitica","Aculeata"))) %>% 
	filter(group != "")

subgroup6 <- egg_database %>% 
	filter(group == "Neuropteroidea") %>% 
	mutate(group = ifelse(order == "Coleoptera",
		suborder,order)) %>% 
	filter(group != "")

subgroup7 <- egg_database %>% 
	filter(group == "Amphiesmenoptera") %>% 
	mutate(group = ifelse(order == "Trichoptera",
		order,ifelse(superfamily %in% c("Papilionoidea","Bombycoidea","Noctuoidea"),
			superfamily,"other moths"))) %>% 
	filter(group != "")

subgroup8 <- egg_database %>% 
	filter(group == "Antliophora") %>% 
	mutate(group = ifelse(order == "Diptera",
		suborder,order)) %>% 
	filter(group != "")

plot_subgroup <- function(subgroup) {
	subgroup_levels <- subgroup %>% group_by(group) %>% summarize(count = n()) %>% arrange(desc(count)) %>% pull(group)
	names(submrk) <- subgroup_levels
	sg1 <- ggplot(subgroup,aes(x = ar, y = vol,color = group)) + 
		geom_hline(yintercept = 10^-2,linetype = "dashed",color = "dark gray",size=1) + 
		geom_vline(xintercept = 2,linetype = "dashed",color = "dark gray",size=1) + 
		geom_point(size=2.25,alpha = 0.75,pch=16) +	
		scale_x_log10(breaks = ar_breaks) + 
		scale_y_log10(breaks = vol_breaks,expand=c(0,0.2)) +
		scale_color_manual(values = submrk) + 
		theme(legend.position = "none")
	return(sg1)
}

# Build theoretical morphospace plots

morpho_aspect_ratio_asymmetry_untransformed <- ggplot(egg_database,aes(x = ar, y = asym,color = group)) + 
	geom_point(alpha = 0.65,shape=16,size=1) + 
	theme(legend.position="none") + 
	scale_color_manual(values = mrk) + 
	scale_x_log10(breaks = ar_breaks) + 
	scale_y_continuous(breaks = asym_breaks_untransformed,expand=c(0,0.2)) + 
	coord_fixed()

morpho_aspect_ratio_curvature_untransformed <- ggplot(egg_database,aes(x = ar, y = curv,color = group)) + 
	geom_point(alpha = 0.65,shape=16,size=1) + 
	theme(legend.position="none") + 
	scale_color_manual(values = mrk) + 
	scale_x_log10(breaks = ar_breaks) + 
	scale_y_continuous(breaks = curv_breaks_untransformed,expand=c(0,0.2)) + 
	coord_fixed(ratio = (1/ 3))

morpho_asymmetry_curvature_untransformed <- ggplot(egg_database,aes(x = asym, y = curv,color = group)) + 
	geom_point(alpha = 0.65,shape=16,size=1) + 
	theme(legend.position="none") + 
	scale_color_manual(values = mrk) + 
	scale_x_continuous(breaks = asym_breaks_untransformed,expand=c(0,0.2)) + 
	scale_y_continuous(breaks = curv_breaks_untransformed,expand=c(0,0.2))  + 
	coord_fixed(ratio=(1/3.49))

# Build morphospaces comparing insect and avian data
insect_bird_asym_ar <- insect_bird_egg_data %>% select(group,sqasym,logar,asym,ar) %>% na.omit()
# build a geometric hull showing the space occupied by the avian egg data
ar_asym_hulls <- plyr::ddply(insect_bird_asym_ar, "group", function(df) df[chull(df$asym, df$ar), ])
ar_vol_hulls <- plyr::ddply(egg_database %>% select(logvol,logar,group) %>% na.omit(), "group", function(df) df[chull(df$logvol, df$logar), ])

morpho_aspect_ratio_asymmetry_birds_untransformed <- ggplot(egg_database,aes(x = asym, y = ar,color = group,fill = group)) + 
	geom_point(alpha = 0.9,pch=16) + 
	theme(legend.position="none") + 
	scale_y_log10(breaks = ar_breaks) + 
	scale_x_continuous(breaks = asym_breaks_untransformed) + 
	scale_color_manual(values = mrk) + 
	geom_polygon(data = ar_asym_hulls %>% filter(group == "Aves"),alpha=0.5,size=1) + 
	scale_fill_manual(values = mrk)


# Print morphospaces
pdf(file="subgroup1.pdf",width=6,height=6,useDingbats =F)
print(plot_subgroup(subgroup1))
dev.off()
pdf(file="subgroup2.pdf",width=6,height=6,useDingbats =F)
print(plot_subgroup(subgroup2))
dev.off()
pdf(file="subgroup3.pdf",width=6,height=6,useDingbats =F)
print(plot_subgroup(subgroup3))
dev.off()
pdf(file="subgroup4.pdf",width=6,height=6,useDingbats =F)
print(plot_subgroup(subgroup4))
dev.off()
pdf(file="subgroup5.pdf",width=6,height=6,useDingbats =F)
print(plot_subgroup(subgroup5))
dev.off()
pdf(file="subgroup6.pdf",width=6,height=6,useDingbats =F)
print(plot_subgroup(subgroup6))
dev.off()
pdf(file="subgroup7.pdf",width=6,height=6,useDingbats =F)
print(plot_subgroup(subgroup7))
dev.off()
pdf(file="subgroup8.pdf",width=6,height=6,useDingbats =F)
print(plot_subgroup(subgroup8))
dev.off()

sink(file="subgroup_legends.txt")
cat("\nApterygota subgroups\n")
print(subgroup1 %>% group_by(group) %>% summarize(count = n()) %>% arrange(desc(count)) %>% pull(group))
cat("\nPalaeoptera subgroups\n")
print(subgroup2 %>% group_by(group) %>% summarize(count = n()) %>% arrange(desc(count)) %>% pull(group))
cat("\nPolyneoptera subgroups\n")
print(subgroup3 %>% group_by(group) %>% summarize(count = n()) %>% arrange(desc(count)) %>% pull(group))
cat("\nCondylognatha subgroups\n")
print(subgroup4 %>% group_by(group) %>% summarize(count = n()) %>% arrange(desc(count)) %>% pull(group))
cat("\nHymenoptera subgroups\n")
print(subgroup5 %>% group_by(group) %>% summarize(count = n()) %>% arrange(desc(count)) %>% pull(group))
cat("\nNeuropteroidea subgroups\n")
print(subgroup6 %>% group_by(group) %>% summarize(count = n()) %>% arrange(desc(count)) %>% pull(group))
cat("\nAmphiesmenoptera subgroups\n")
print(subgroup7 %>% group_by(group) %>% summarize(count = n()) %>% arrange(desc(count)) %>% pull(group))
cat("\nAntliophora subgroups\n")
print(subgroup8 %>% group_by(group) %>% summarize(count = n()) %>% arrange(desc(count)) %>% pull(group))
sink()


pdf(file="polyembryonic.pdf",width=6,height=6,useDingbats =F)
print(polyembryonic_plot)
dev.off()

pdf(file="ootheca.pdf",width=6,height=6,useDingbats =F)
print(ootheca_plot)
dev.off()

pdf(file="morpho_aspect_ratio_asymmetry.pdf",width=6,height=6,useDingbats =F)
print(morpho_aspect_ratio_asymmetry_untransformed)
dev.off()

pdf(file="morpho_aspect_ratio_curvature.pdf",width=6,height=6,useDingbats =F)
print(morpho_aspect_ratio_curvature_untransformed)
dev.off()

pdf(file="morpho_asymmetry_curvature.pdf",width=6,height=6,useDingbats =F)
print(morpho_asymmetry_curvature_untransformed)
dev.off()

pdf(file="morpho_aspect_ratio_asymmetry_birds.pdf",width=6,height=6,useDingbats =F)
print(morpho_aspect_ratio_asymmetry_birds_untransformed)
dev.off()

# Build strip plots comparing occupancy of parameter space by insect and avian groups
# Set the order of groups in strip plots
insect_bird_egg_data$group <- factor(insect_bird_egg_data$group, levels = rev(group_levels))

vol_strip <- ggplot(insect_bird_egg_data,aes(x = group, y = vol, color = group)) + 
	geom_jitter(position = position_jitter(0.45)) + 
	scale_color_manual(values = mrk) + 
	scale_y_log10(breaks = vol_breaks) + 
	coord_flip() + 
	theme(legend.position = "none")

ar_strip <- ggplot(insect_bird_egg_data,aes(x = group, y = ar, color = group)) + 
	geom_jitter(position = position_jitter(0.45)) + 
	scale_color_manual(values = mrk) + 
	scale_y_log10(breaks = ar_breaks) + 
	coord_flip() + 
	theme(legend.position = "none")

asym_strip_untransformed <- ggplot(insect_bird_egg_data,aes(x = group, y = asym, color = group)) + 
	geom_jitter(position = position_jitter(0.45)) + 
	scale_color_manual(values = mrk) + 
	scale_y_continuous(breaks = asym_breaks_untransformed) + 
	coord_flip() + 
	theme(legend.position = "none")

curv_strip_untransformed <- ggplot(insect_bird_egg_data  %>% filter(!(group == "Aves")),aes(x = group, y = curv, color = group)) + 
	geom_jitter(position = position_jitter(0.45)) + 
	scale_color_manual(values = mrk) + 
	scale_y_continuous(breaks = curv_breaks_untransformed) + 
	coord_flip() + 
	theme(legend.position = "none")

# Print strip plots
pdf(file="morpho_volume_aspect_ratio.pdf",width=9.3,height=8,useDingbats =F)
print(morpho_volume_aspect_ratio)
dev.off()
pdf(file="ar_strip.pdf",width=12,height=8,useDingbats =F)
print(ar_strip)
dev.off()
pdf(file="vol_strip.pdf",width=12,height=8,useDingbats=F)
print(vol_strip)
dev.off()
pdf(file="curv_strip.pdf",width=12,height=8,useDingbats =F)
print(curv_strip_untransformed)
dev.off()
pdf(file="asym_strip.pdf",width=12,height=8,useDingbats =F)
print(asym_strip_untransformed)
dev.off()
