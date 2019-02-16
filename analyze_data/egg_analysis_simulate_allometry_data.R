### Created by SHC April 4, 2018

### This code simulates a new egg database with a given allometric relationship

source("analyze_data/egg_analysis_build_dataframe.R")
source("analyze_data/egg_analysis_read_trees.R")
source("analyze_data/egg_analysis_pgls_functions.R")
source("analyze_data/egg_analysis_geiger_functions.R")
library(phylolm)
library(phytools)

### Simulate data with Phylolm

# The slope of regression between length and width for simulation is set by an argument
args = commandArgs(trailingOnly=TRUE)
slope <- as.numeric(args[1])

# Pick a tree for model fitting (PGLS will be run over posterior distribution)
tree <- genus_mcc_tree

# Fit an evolutionary model to the data
# trait1 = egg length
ms_X1 <- egg_database %>% mutate(trait = logX1,rank = genus) %>% build_geiger_dataset(tree = tree)
ms_X1 <- ms_X1[order(match(names(ms_X1),tree$tip.label))]
ms_X1_trim <- drop.tip(tree,setdiff(tree$tip.label,names(ms_X1)))
X1_fit <- fitContinuous(ms_X1_trim,ms_X1,model="EB")

# trait2 = egg width
ms_X2 <- egg_database %>% mutate(trait = logX2, rank = genus) %>% build_geiger_dataset(tree = tree)
ms_X2 <- ms_X2[order(match(names(ms_X2),tree$tip.label))]
ms_X2_trim <- drop.tip(tree,setdiff(tree$tip.label,names(ms_X2)))
X2_fit <- fitContinuous(ms_X2_trim,ms_X2,model="EB")

# set groups for subgroup analyses
group_list <- c("Hymenoptera","Condylognatha","Polyneoptera","Palaeoptera","Amphiesmenoptera","Antliophora","Neuropteroidea") 

### This function uses the fitted model parameters to simulate new trait data
### It performs a PGLS comparison of the new data and can be repeated
run_all_allometry_pgls <- function(tree){
	# simulate new length based on fitted parameters
	x_rate <- X1_fit$opt$sigsq
	x_anc <- X1_fit$opt$z0
	x_a <- X1_fit$opt$a
	x <- x_anc + rTrait(n = 1,phy = unroot(tree),model = "EB",parameters = list(a = x_a, sigma2 = x_rate))

	# simulate new width with known slope based on fitted parameters
	y_rate <- X2_fit$opt$sigsq
	y_anc <- X2_fit$opt$z0
	y_a<- X2_fit$opt$a
	# slope is set by slope*x
	y <- y_anc + slope*x + rTrait(n = 1,phy = unroot(tree),model = "EB",parameters = list(a = y_a, sigma2 = y_rate))
	# slope of 0, uncomment here
	#y <- y_anc + rTrait(n = 1,phy = unroot(tree),model = "EB",parameters = list(a = y_a, sigma2 = y_rate))
	 
	# combine the simulated datasets
	data_of_interest <- data.frame(genus = tree$tip.label, simlogX1 = x, simlogX2 = y)
	egg_database_sim <- merge(egg_database,data_of_interest,by = "genus")

	# remove taxa with NA data in database - mask missing data
	X1_dat <- data.frame(trait1 = egg_database$X1, rank = egg_database$genus) %>% na.omit()
	X2_dat <- data.frame(trait1 = egg_database$X2, rank = egg_database$genus) %>% na.omit()
	egg_database_sim <- egg_database_sim %>% filter(genus %in% X2_dat$rank) %>% filter(genus %in% X1_dat$rank)

	# Run a PGLS on simulated length vs width
	length_width <- egg_database_sim %>% select(genus,simlogX1,simlogX2,group) %>% rename(rank = genus, trait1 = simlogX1, trait2= simlogX2)
	allometry_pgls_length_width <- run_all_taxa_and_by_group_pgls(length_width,tree,"sim_length_width",group_list)

	return(list(allometry_pgls_length_width))
}

### Repeat the simulation over the posterior distribution of trees
allometry_distribution_raw <- lapply(genus_trees,run_all_allometry_pgls)

### Create a table of results across all insects
allometry_all_taxa_table <- round(data.frame(
	slope_min = c(
	min(sapply(allometry_distribution_raw,function(x) {x[[1]][[1]]$coefficients[[2]]}))),
	slope_max = c(
	max(sapply(allometry_distribution_raw,function(x) {x[[1]][[1]]$coefficients[[2]]}))),
	int_min = c(
	min(sapply(allometry_distribution_raw,function(x) {x[[1]][[1]]$coefficients[[1]]}))),
	int_max = c(
	max(sapply(allometry_distribution_raw,function(x) {x[[1]][[1]]$coefficients[[1]]}))),
	pval_min = c(
	min(sapply(allometry_distribution_raw,function(x) {x[[1]][[1]]$tTable[[8]]}))),
	pval_max = c(
	max(sapply(allometry_distribution_raw,function(x) {x[[1]][[1]]$tTable[[8]]}))),
	sample_size = c(
	min(sapply(allometry_distribution_raw,function(x) {x[[1]][[1]]$dims$N}))),
	row.names = c("length vs width")),
	4)

### Create a table of results by group
get_group_allometry_table <- function(value) {
	table <- round(data.frame(
	slope_min = c(
	min(sapply(allometry_distribution_raw,function(x) {x[[1]][[2]][[value]]$coefficients[[2]]}))),
	slope_max = c(
	max(sapply(allometry_distribution_raw,function(x) {x[[1]][[2]][[value]]$coefficients[[2]]}))),
	int_min = c(
	min(sapply(allometry_distribution_raw,function(x) {x[[1]][[2]][[value]]$coefficients[[1]]}))),
	int_max = c(
	max(sapply(allometry_distribution_raw,function(x) {x[[1]][[2]][[value]]$coefficients[[1]]}))),
	pval_min = c(
	min(sapply(allometry_distribution_raw,function(x) {x[[1]][[2]][[value]]$tTable[[8]]}))),
	pval_max = c(
	max(sapply(allometry_distribution_raw,function(x) {x[[1]][[2]][[value]]$tTable[[8]]}))),
	sample_size = c(
	min(sapply(allometry_distribution_raw,function(x) {x[[1]][[2]][[value]]$dims$N}))),
	row.names = c("length vs width")),
	4)
	return(table)
}

### Compare the slope of allometric regressions across phylogenetic clades
# This function gets the slopes from a distribution of regression results
get_slope_distribution <- function(pgls,value) {
	slope <- sapply(allometry_distribution_raw,function(x) {x[[1]][[2]][[value]]$coefficients[[2]]})
	return(slope)
}

# Get the slopes over the distribution of results
slope_dist_l_w <- lapply(seq(length(group_list)),get_slope_distribution,pgls=1)
names(slope_dist_l_w) <- group_list
slope_dist_l_w <- reshape::melt(slope_dist_l_w)

### Get results table for each group
allometry_Hymenoptera_table <- get_group_allometry_table(which(group_list == "Hymenoptera"))
allometry_Condylognatha_table <- get_group_allometry_table(which(group_list == "Condylognatha"))
allometry_Neuropteroidea_table <- get_group_allometry_table(which(group_list == "Neuropteroidea"))
allometry_Amphiesmenoptera_table <- get_group_allometry_table(which(group_list == "Amphiesmenoptera"))
allometry_Antliophora_table <- get_group_allometry_table(which(group_list == "Antliophora"))
allometry_Palaeoptera_table <- get_group_allometry_table(which(group_list == "Palaeoptera"))
allometry_Polyneoptera_table <- get_group_allometry_table(which(group_list == "Polyneoptera"))

### Print the results
sink(file=paste("slope",slope,"sim_allometry_tables.txt",sep="_"))
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

save.image(file=paste("slope",slope,"egg_analysis_simulate_allometry_workspace.RData",sep="_"))
