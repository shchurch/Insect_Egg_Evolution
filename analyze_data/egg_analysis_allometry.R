### Created by SHC April 2018

### This code executes the allometry comparisons on the insect egg data

source("analyze_data/egg_analysis_build_dataframe.R")
source("analyze_data/egg_analysis_read_trees.R")
library(phytools)

analysis_name <- ""
args = commandArgs(trailingOnly=TRUE)
# This program takes two arguments
# 1 - which backbone tree to use, 'misof' [default] or 'rainford'
# 2 - which correlation matrix to use 'brownian' [default] or 'blomberg'
if(args[1] == "rainford") {
	genus_trees <- rainford_genus_trees
	genus_mcc_tree <- rainford_genus_mcc_tree
	analysis_name <- paste("rainford",analysis_name,sep="_")
}
if(args[2] == "blomberg") {
	blom_flag <- TRUE
	analysis_name <- paste("corBlomberg",analysis_name,sep="_")
}
source("analyze_data/egg_analysis_pgls_functions.R")

# Set factor for resampling the body size dataframe
downsample_factor <- 0.5
fam_count_threshold <- 1

### set groups for subgroup analyses
group_list <- c("Hymenoptera","Condylognatha","Antliophora","Neuropteroidea","Amphiesmenoptera","Polyneoptera","Palaeoptera") 

### This function performs 5 PGLS comparisons of 2 traits, and returns the results as a list of 2 lists
### Each analysis is performed over all insects (list 1), and across each group (list 2)
### This function can be run over a distribution of trees
run_all_allometry_pgls <- function(tree){
	### PGLS length vs width
	length_width <- egg_database %>% select(genus,logX1,logX2,group) %>% rename(rank = genus, trait1 = logX1, trait2= logX2)
	allometry_pgls_length_width <- run_all_taxa_and_by_group_pgls(length_width,tree,paste(analysis_name,"length_width",sep=""),group_list)

	### PGLS length vs asymmetry, width = independent
	length_asym_width <- egg_database %>% select(genus,logX1,sqasym,logX2,group) %>% rename(rank = genus, trait1 = logX1, trait2= sqasym, indep = logX2)
	allometry_pgls_length_asym_width <- run_all_taxa_and_by_group_resid_pgls(length_asym_width,tree,paste(analysis_name,"length_asym_width",sep=""),group_list)

	### PGLS length vs curvature, width = independent
	length_curv_width <- egg_database %>% select(genus,logX1,sqcurv,logX2,group) %>% rename(rank = genus, trait1 = logX1, trait2= sqcurv, indep = logX2)
	allometry_pgls_length_curv_width <- run_all_taxa_and_by_group_resid_pgls(length_curv_width,tree,paste(analysis_name,"length_curv_width",sep=""),group_list)

	### build body size dataset
	source("analyze_data/egg_analysis_body_size.R")

	### PGLS egg size vs body size
	vol_body <- egg_database_family_body %>% select(family,logvol,logbodyvol,group) %>% rename(rank = family, trait1 = logbodyvol, trait2= logvol)
	allometry_pgls_vol_body <- run_all_taxa_and_by_group_pgls(vol_body,fam_tree,paste(analysis_name,"vol_body",sep=""),group_list)

	### PGLS length vs width, body size = independent
	length_width_body <- egg_database_family_body %>% select(family,logX1,logX2,logbody,group) %>% rename(rank = family, trait1 = logX1, trait2= logX2, indep = logbody)
	allometry_pgls_length_width_body <- run_all_taxa_and_by_group_resid_pgls(length_width_body,fam_tree,paste(analysis_name,"length_width_body",sep=""),group_list)

	return(list(allometry_pgls_length_width,allometry_pgls_length_asym_width,allometry_pgls_length_curv_width,allometry_pgls_vol_body,allometry_pgls_length_width_body))
}

### Run all PGLS comparisons over the posterior distribution of trees
### Results stored as a list of lists
allometry_distribution_raw <- lapply(genus_trees,run_all_allometry_pgls)

### Build a results table for all insects
get_allometry_all_taxa_table <- function(pgls) {
	table <- round(data.frame(
	slope_min = min(sapply(allometry_distribution_raw,function(x) {x[[pgls]][[1]]$coefficients[[2]]})),
	slope_max = max(sapply(allometry_distribution_raw,function(x) {x[[pgls]][[1]]$coefficients[[2]]})),
	int_min = min(sapply(allometry_distribution_raw,function(x) {x[[pgls]][[1]]$coefficients[[1]]})),
	int_max = max(sapply(allometry_distribution_raw,function(x) {x[[pgls]][[1]]$coefficients[[1]]})),
	pval_min = min(sapply(allometry_distribution_raw,function(x) {x[[pgls]][[1]]$tTable[[8]]})),
	pval_max = max(sapply(allometry_distribution_raw,function(x) {x[[pgls]][[1]]$tTable[[8]]})),
	sample_size = min(sapply(allometry_distribution_raw,function(x) {x[[pgls]][[1]]$dims$N}))),
	4)
	return(table)
}

allometry_l_w_table <- get_allometry_all_taxa_table(1)
allometry_l_asym_w_table <- get_allometry_all_taxa_table(2)
allometry_l_curv_w_table <- get_allometry_all_taxa_table(3)
allometry_vol_body_table <- get_allometry_all_taxa_table(4)
allometry_l_w_body_table <- get_allometry_all_taxa_table(5)

### Build a results table for each group
get_group_allometry_table <- function(pgls,value) {
	table <- round(data.frame(
	slope_min = min(sapply(allometry_distribution_raw,function(x) {x[[pgls]][[2]][[value]]$coefficients[[2]]})),
	slope_max = max(sapply(allometry_distribution_raw,function(x) {x[[pgls]][[2]][[value]]$coefficients[[2]]})),
	int_min = min(sapply(allometry_distribution_raw,function(x) {x[[pgls]][[2]][[value]]$coefficients[[1]]})),
	int_max = max(sapply(allometry_distribution_raw,function(x) {x[[pgls]][[2]][[value]]$coefficients[[1]]})),
	pval_min = min(sapply(allometry_distribution_raw,function(x) {x[[pgls]][[2]][[value]]$tTable[[8]]})),
	pval_max = max(sapply(allometry_distribution_raw,function(x) {x[[pgls]][[2]][[value]]$tTable[[8]]})),
	sample_size = min(sapply(allometry_distribution_raw,function(x) {x[[pgls]][[2]][[value]]$dims$N}))),
	4)
	return(table)
}

### Get the results table for each group
allometry_group_l_w_table <- lapply(seq(length(group_list)),get_group_allometry_table,pgls=1)
allometry_group_l_asym_w_table <- lapply(seq(length(group_list)),get_group_allometry_table,pgls=2)
allometry_group_l_curv_w_table <- lapply(seq(length(group_list)),get_group_allometry_table,pgls=3)
allometry_group_vol_body_table <- lapply(seq(length(group_list)),get_group_allometry_table,pgls=4)
allometry_group_l_w_body_table <- lapply(seq(length(group_list)),get_group_allometry_table,pgls=5)

### Compare the slope of allometric regressions across phylogenetic clades
# This function gets the slopes from a distribution of regression results
get_slope_distribution <- function(pgls,value) {
	slope <- sapply(allometry_distribution_raw,function(x) {x[[pgls]][[2]][[value]]$coefficients[[2]]})
	return(slope)
}

# Get the slopes over the distribution of results
slope_dist_l_w <- lapply(seq(length(group_list)),get_slope_distribution,pgls=1)
# name them according to the ecological regimes
names(slope_dist_l_w) <- group_list
# reshape the list for comparison
slope_dist_l_w <- reshape::melt(slope_dist_l_w)
# plot the comparison
slope_l_w_plot <- ggplot(slope_dist_l_w,aes(x = factor(L1,levels=group_levels), y = value, color = L1, fill = L1)) + 
	geom_abline(intercept=1.0,slope=0,linetype=2) + 
	geom_boxplot(outlier.colour = NULL) + 
	scale_fill_manual(values = mrk) + 
	scale_color_manual(values = mrk) + 
	stat_summary(geom = "crossbar", width=0.65, lwd = 0.7, fatten=0, color="white", fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x))) }) + 	
	theme(legend.position="none",axis.text.x = element_blank()) + ylab("slope") + xlab("clade") + ylim(c(0.5,1.0))

pdf(paste(analysis_name,"l_w_allometry_slope_plot.pdf",sep=""),height=3.5,width=7)
print(slope_l_w_plot)
dev.off()

# Get the slopes over the distribution of results
slope_dist_vol_body_vol <- lapply(seq(length(group_list)),get_slope_distribution,pgls=4)
# name them according to the ecological regimes
names(slope_dist_vol_body_vol) <- group_list
# reshape the list for comparison
slope_dist_vol_body_vol <- reshape::melt(slope_dist_vol_body_vol)
# plot the comparison
slope_vol_body_vol_plot <- ggplot(slope_dist_vol_body_vol,aes(x = factor(L1,levels=group_levels), y = value, color = L1, fill = L1)) + 
	geom_abline(intercept=1.0,slope=0,linetype=2) + 
	geom_boxplot(outlier.colour = NULL) + 
	scale_fill_manual(values = mrk) + 
	scale_color_manual(values = mrk) + 
	stat_summary(geom = "crossbar", width=0.65, lwd = 0.7, fatten=0, color="white", fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x))) }) + 	
	theme(legend.position="none",axis.text.x = element_blank()) + ylab("slope") + xlab("clade") + ylim(c(0.2,1.0))

pdf(paste(analysis_name,"vol_body_vol_allometry_slope_plot.pdf",sep=""),height=3.5,width=7)
print(slope_vol_body_vol_plot)
dev.off()

# Get the slopes over the distribution of results
slope_dist_l_w_body <- lapply(seq(length(group_list)),get_slope_distribution,pgls=5)
# name them according to the ecological regimes
names(slope_dist_l_w_body) <- group_list
# reshape the list for comparison
slope_dist_l_w_body <- reshape::melt(slope_dist_l_w)
# plot the comparison
slope_l_w_body_plot <- ggplot(slope_dist_l_w_body,aes(x = factor(L1,levels=group_levels), y = value, color = L1, fill = L1)) + 
	geom_abline(intercept=1.0,slope=0,linetype=2) + 
	geom_boxplot(outlier.colour = NULL) + 
	scale_fill_manual(values = mrk) + 
	scale_color_manual(values = mrk) + 
	stat_summary(geom = "crossbar", width=0.65, lwd = 0.7, fatten=0, color="white", fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x))) }) + 	
	theme(legend.position="none",axis.text.x = element_blank()) + ylab("slope") + xlab("clade") + ylim(c(0.6,1.0))

pdf(paste(analysis_name,"l_w_body_allometry_slope_plot.pdf",sep=""),height=4,width=7)
print(slope_l_w_body_plot)
dev.off()

### Print the results
sink(paste(analysis_name,"allometry_tables.txt",sep=""))
cat("Across all insects\n")
cat("length_width")
print(allometry_l_w_table)
cat("length_asymmetry_width")
print(allometry_l_asym_w_table)
cat("length_curvature_width")
print(allometry_l_curv_w_table)
cat("egg_volume_body")
print(allometry_vol_body_table)
cat("length_width_body")
print(allometry_l_w_body_table)
for(i in seq(1:length(group_list))) {
	cat("\n")
	print(group_list[i])
	cat("\nlength_width")
	print(allometry_group_l_w_table[[i]])
	cat("length_asymmetry_width")
	print(allometry_group_l_asym_w_table[[i]])
	cat("length_curvature_width")
	print(allometry_group_l_curv_w_table[[i]])
	cat("egg_volume_body")
	print(allometry_group_vol_body_table[[i]])
	cat("length_width_body")
	print(allometry_group_l_w_body_table[[i]])
}

save.image(file=paste("egg_analysis_",analysis_name,"allometry_workspace.RData",sep=""))


