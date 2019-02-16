### Created by SHC April 2018

### This code simulates new egg + body size allometry data with a fixed slope

source("analyze_data/egg_analysis_build_dataframe.R")
# Set factor for resampling the body size dataframe, here no downsampling used
downsample_factor <- 1.0
fam_count_threshold <- 0
source("analyze_data/egg_analysis_body_size.R")
source("analyze_data/egg_analysis_read_trees.R")
source("analyze_data/egg_analysis_pgls_functions.R")
source("analyze_data/egg_analysis_geiger_functions.R")
library(phylolm)
library(phytools)

### Simulate data with Phylolm

# The slope of regression between length and width for simulation is set by an argument
args = commandArgs(trailingOnly=TRUE)
slope <- as.numeric(args[1])

# Use the family level tree
tree <- fam_tree

# Fit an evolutionary model to the data
# trait1 = body size
ms_body <- egg_database_family_body %>% mutate(trait = logbodyvol,rank = family) %>% build_geiger_dataset(tree = tree)
ms_body <- ms_body[order(match(names(ms_body),tree$tip.label))]
ms_body_trim <- drop.tip(tree,setdiff(tree$tip.label,names(ms_body)))
body_fit <- fitContinuous(ms_body_trim,ms_body,model="OU")

# trait2 = egg width
ms_vol <- egg_database_family_body %>% mutate(trait = logvol, rank = family) %>% build_geiger_dataset(tree = tree)
ms_vol <- ms_vol[order(match(names(ms_vol),tree$tip.label))]
ms_vol_trim <- drop.tip(tree,setdiff(tree$tip.label,names(ms_vol)))
vol_fit <- fitContinuous(ms_vol_trim,ms_vol,model="BM")

### set groups for subgroup analyses
group_list <- c("Hymenoptera","Condylognatha","Polyneoptera","Palaeoptera","Amphiesmenoptera","Antliophora","Neuropteroidea") 

### This function uses the fitted model parameters to simulate new trait data
### It performs a PGLS comparison of the new data and can be repeated
run_body_allometry_pgls <- function(tree){
	# simulate new body size based on fitted parameters
	b_rate <- body_fit$opt$sigsq
	b_anc <- body_fit$opt$z0
	b_alpha <- body_fit$opt$alpha
	b <- b_anc + rTrait(n = 1,phy = unroot(tree),model = "OU",parameters = list(alpha = b_alpha, sigma2 = b_rate, ancestral.state = b_anc))

	# simulate new egg volume with known slope based on fitted parameters
	v_rate <- vol_fit$opt$sigsq
	v_anc <- vol_fit$opt$z0
	# slope is set by slope*b
	v <- v_anc + slope*b + rTrait(n = 1,phy = unroot(tree),model = "BM",parameters = list(ancestral.state = v_anc, sigma2 = v_rate))
	# slope of 0, uncomment here
	#v <- v_anc + rTrait(n = 1,phy = unroot(tree),model = "BM",parameters = list(ancestral.state = v_anc, sigma2 = v_rate))
	 
	# combine the simulated datasets
	data_of_interest <- data.frame(family = tree$tip.label, simlogbodyvol = b, simlogvol = v)
	egg_database_family_body_sim <- merge(egg_database_family_body,data_of_interest,by = "family")

	# remove taxa with NA data in database - mask missing data
	body_dat <- data.frame(trait1 = egg_database_family_body$logbodyvol, rank = egg_database_family_body$family) %>% na.omit()
	vol_dat <- data.frame(trait1 = egg_database_family_body$logvol, rank = egg_database_family_body$family) %>% na.omit()
	egg_database_family_body_sim <- egg_database_family_body_sim %>% filter(family %in% vol_dat$rank) %>% filter(family %in% body_dat$rank)

	# Run a PGLS on simulated body size vs egg size
	bodyvol_vol <- egg_database_family_body_sim %>% select(family,simlogbodyvol,simlogvol,group) %>% rename(rank = family, trait1 = simlogbodyvol, trait2= simlogvol)
	allometry_pgls_bodyvol_vol <- run_all_taxa_and_by_group_pgls(bodyvol_vol,tree,"sim_bodyvol_vol",group_list)

	return(list(allometry_pgls_bodyvol_vol))
}

### Repeat the simulation 100 times, simulating new data each time
allometry_distribution_raw <- replicate(100,run_body_allometry_pgls(fam_tree))

### Create a table of results across all insects
allometry_all_taxa_table <- round(data.frame(
	slope_min = c(
	min(sapply(allometry_distribution_raw,function(x) {x[[1]]$coefficients[[2]]}))),
	slope_max = c(
	max(sapply(allometry_distribution_raw,function(x) {x[[1]]$coefficients[[2]]}))),
	pval_min = c(
	min(sapply(allometry_distribution_raw,function(x) {x[[1]]$tTable[[8]]}))),
	pval_max = c(
	max(sapply(allometry_distribution_raw,function(x) {x[[1]]$tTable[[8]]}))),
	sample_size = c(
	min(sapply(allometry_distribution_raw,function(x) {x[[1]]$dims$N}))),
	row.names = c("bodyvol vs eggvol")),
	4)

### Create a table of results by group
get_group_allometry_table <- function(value) {
	table <- round(data.frame(
	slope_min = c(
	min(sapply(allometry_distribution_raw,function(x) {x[[2]][[value]]$coefficients[[2]]}))),
	slope_max = c(
	max(sapply(allometry_distribution_raw,function(x) {x[[2]][[value]]$coefficients[[2]]}))),
	pval_min = c(
	min(sapply(allometry_distribution_raw,function(x) {x[[2]][[value]]$tTable[[8]]}))),
	pval_max = c(
	max(sapply(allometry_distribution_raw,function(x) {x[[2]][[value]]$tTable[[8]]}))),
	sample_size = c(
	min(sapply(allometry_distribution_raw,function(x) {x[[2]][[value]]$dims$N}))),
	row.names = c("bodyvol vs eggvol")),
	4)
	return(table)
}

### Get results table for each group
allometry_Hymenoptera_table <- get_group_allometry_table(which(group_list == "Hymenoptera"))
allometry_Condylognatha_table <- get_group_allometry_table(which(group_list == "Condylognatha"))
allometry_Neuropteroidea_table <- get_group_allometry_table(which(group_list == "Neuropteroidea"))
allometry_Amphiesmenoptera_table <- get_group_allometry_table(which(group_list == "Amphiesmenoptera"))
allometry_Antliophora_table <- get_group_allometry_table(which(group_list == "Antliophora"))
allometry_Palaeoptera_table <- get_group_allometry_table(which(group_list == "Palaeoptera"))
allometry_Polyneoptera_table <- get_group_allometry_table(which(group_list == "Polyneoptera"))

### Compare the slope of allometric regressions across phylogenetic clades
# This function gets the slopes from a distribution of regression results
get_slope_distribution <- function(pgls,value) {
	slope <- sapply(allometry_distribution_raw,function(x) {x[[2]][[value]]$coefficients[[2]]})
	return(slope)
}

# Get the slopes over the distribution of results
slope_dist_vol_bodyvol <- lapply(seq(length(group_list)),get_slope_distribution,pgls=1)
names(slope_dist_vol_bodyvol) <- group_list
slope_dist_vol_bodyvol <- reshape::melt(slope_dist_vol_bodyvol)

### Print the results
sink(paste("slope",slope,"sim_body_tables.txt",sep="_"))
cat("Across all insects")
print(allometry_all_taxa_table)
cat("\nHymenoptera")
print(allometry_Hymenoptera_table)
cat("\nCondylognatha")
print(allometry_Condylognatha_table)
cat("\nNeuropteroidea")
print(allometry_Neuropteroidea_table)
cat("\nAmphiesmenoptera")
print(allometry_Amphiesmenoptera_table)
cat("\nAntliophora")
print(allometry_Antliophora_table)
cat("\nPalaeoptera")
print(allometry_Palaeoptera_table)
cat("\nPolyneoptera")
print(allometry_Polyneoptera_table)
sink()

save.image(file=paste("slope",slope,"egg_analysis_simulate_body_workspace.RData",sep="_"))

