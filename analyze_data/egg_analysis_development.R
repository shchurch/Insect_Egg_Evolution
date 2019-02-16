### Created by SHC April 2018

### This code compares egg size and developmental rate using phylogenetic comparative methods

source("analyze_data/egg_analysis_build_dataframe.R")
source("analyze_data/egg_analysis_read_trees.R")
source("analyze_data/egg_analysis_geiger_functions.R")

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

### Read in the developmental tables
development <- read.delim("analyze_data/development.csv",header=T,stringsAsFactors= F)

### Set up the developmental parameters
development <- development %>% 
		# correct for temperature
		mutate(Tm_hatch = temp_to_hatch + 273.15,
		Tm_development = temp_development + 273.15,
		int_btw_preblast_mitoses_hrs = int_btw_preblast_mitoses / 60,
		Da_int_btw_preblast_mitoses = (int_btw_preblast_mitoses_hrs) / exp((8000/293.15)-(8000/Tm_development)),
		Da_midcellularization = (time_to_midcellularization) / exp((8000/293.15)-(8000/Tm_development)),
		Da_hatch = (time_to_hatch_hours) / exp((8000/293.15)-(8000/Tm_hatch)),
		# log transform the values
		logDa_int_btw_preblast_mitoses = log10(Da_int_btw_preblast_mitoses),
		logDa_midcellularization = log10(Da_midcellularization),
		logDa_hatch = log10(Da_hatch))

### This function combines the development and egg datasets
### maximizing the number of taxa with data in the final dataframe
combine_dev_egg_datasets <- function(egg,dev) {
	# set up the development database
	dev_db <- dev %>% 
		# randomly shuffle the data
		sample_frac(1L) %>% 
		# data are combined by exact species match
		mutate(combo_rank = name) %>% 
		# only allow entries with a recorded species name
		filter(!(species == "")) %>% 
		# select only the traits of interest
		select(combo_rank,trait2,genus) %>% 
		# remove all missing data
		na.omit()

	# set up the egg database
	egg_db <- egg %>% 
		# randomly shuffle the data
		sample_frac(1L) %>% 
		# combine by species name
		mutate(combo_rank = name) %>% 
		# select only the traits of interest
		select(combo_rank,group,trait1) %>% 
		# remove all missing data
		na.omit() %>% 
		# remove all taxa without development data
		filter(combo_rank %in% dev_db$combo_rank) %>%
		# group by the rank of interest 
		group_by(combo_rank) %>% 
		# choose one random representative per rank
		slice(1L)

	# set up the combined dataset
	egg_dev <- dev_db %>%
		# remove all development taxa without egg data
		filter(combo_rank %in% egg_db$combo_rank) %>% 
		# sort and group by the rank of interest
		group_by(combo_rank) %>% 
		# choose one random representative per rank
		slice(1L) %>% 
		# join together the development and egg databases
		left_join(egg_db,by="combo_rank")

	return(egg_dev)
}

### Run the PGLS comparisons of developmental traits of interest
run_egg_dev_pgls <- function(tree) {
	# compare egg volume and time to hatching - temperature corrected
	egg_vol <- egg_database %>% mutate(trait1 = logvol)
	dev_hatch <- development %>% mutate(trait2 = logDa_hatch) %>% filter(dipause != "y")
	dat_vol_Da_hatch <- combine_dev_egg_datasets(egg_vol,dev_hatch) %>% rename(rank = genus) %>% build_pgls(tree = tree)
	vol_Da_hatch_pgls_out <- pgls.trees(dat_vol_Da_hatch,genus_mcc_tree,"trait2 ~ trait1")

	# compare egg volume and the interval between mitoses - temperature corrected
	egg_vol <- egg_database %>% mutate(trait1 = logvol)
	dev_Da_int_btw_preblast_mitoses <- development %>% mutate(trait2 = logDa_int_btw_preblast_mitoses) %>% filter(dipause != "y")
	dat_vol_Da_int_btw_preblast_mitoses <- combine_dev_egg_datasets(egg_vol,dev_Da_int_btw_preblast_mitoses) %>% rename(rank = genus) %>% build_pgls(tree = tree)
	vol_Da_int_btw_preblast_mitoses <- pgls.trees(dat_vol_Da_int_btw_preblast_mitoses,genus_mcc_tree,"trait2 ~ trait1")

	# compare egg volume and time to midcellularization  - temperature corrected
	egg_vol <- egg_database %>% mutate(trait1 = logvol)
	dev_Da_time_midcellularization <- development %>% mutate(trait2 = logDa_midcellularization) %>% filter(dipause != "y")
	dat_vol_Da_time_midcellularization <- combine_dev_egg_datasets(egg_vol,dev_Da_time_midcellularization) %>% rename(rank = genus) %>% build_pgls(tree = tree)
	vol_Da_time_midcellularization <- pgls.trees(dat_vol_Da_time_midcellularization,genus_mcc_tree,"trait2 ~ trait1")

	return(list(vol_Da_hatch_pgls_out,vol_Da_int_btw_preblast_mitoses,vol_Da_time_midcellularization))
}

dev_allometry_distribution_raw <- lapply(genus_trees,run_egg_dev_pgls)

### Build a results table for all insects
get_allometry_all_taxa_table <- function(pgls) {
	table <- round(data.frame(
	slope_min = min(sapply(dev_allometry_distribution_raw,function(x) {x[[pgls]]$coefficients[[2]]})),
	slope_max = max(sapply(dev_allometry_distribution_raw,function(x) {x[[pgls]]$coefficients[[2]]})),
	int_min = min(sapply(dev_allometry_distribution_raw,function(x) {x[[pgls]]$coefficients[[1]]})),
	int_max = max(sapply(dev_allometry_distribution_raw,function(x) {x[[pgls]]$coefficients[[1]]})),
	pval_min = min(sapply(dev_allometry_distribution_raw,function(x) {x[[pgls]]$tTable[[8]]})),
	pval_max = max(sapply(dev_allometry_distribution_raw,function(x) {x[[pgls]]$tTable[[8]]})),
	sample_size = min(sapply(dev_allometry_distribution_raw,function(x) {x[[pgls]]$dims$N}))),
	4)
	return(table)
}

# Get the results table by PGLS comparison
allometry_vol_Da_hatch_table <- get_allometry_all_taxa_table(1)
allometry_vol_Da_int_btw_preblast_mitoses_table <- get_allometry_all_taxa_table(2)
allometry_vol_Da_time_midcellularization <- get_allometry_all_taxa_table(3)

### Print the results
sink(paste(analysis_name,"development_tables.txt",sep=""))
cat("Across all insects\n")
cat("volume_time_to_hatching")
print(allometry_vol_Da_hatch_table)
cat("volume_interval_btw_preblast_mitoses")
print(allometry_vol_Da_int_btw_preblast_mitoses_table)
cat("volume_time_midcellularization")
print(allometry_vol_Da_time_midcellularization)
sink()

save.image(file=paste("egg_analysis_",analysis_name,"development_workspace.RData",sep=""))

### Plot the comparisons

# compare egg volume and time to hatching - temperature corrected
egg_vol <- egg_database %>% mutate(trait1 = logvol)
dev_hatch <- development %>% mutate(trait2 = logDa_hatch) %>% filter(dipause != "y")
dat_vol_Da_hatch <- combine_dev_egg_datasets(egg_vol,dev_hatch)
dat_vol_Da_hatch_plot <- ggplot(dat_vol_Da_hatch,aes(x = trait1, y = trait2, color = group)) + 
						geom_point() + 
						scale_color_manual(values = mrk) + 
						theme(legend.position = "none")

pdf(paste(analysis_name,"dat_vol_Da_hatch.pdf",sep=""),width=4,height=4,useDingbats=F)
print(dat_vol_Da_hatch_plot)
dev.off()

# compare egg volume and time to hatching with no phylogenetic correction
non_phylo_regression <- ggplot(dat_vol_Da_hatch,aes(x = trait1, y = trait2)) + 
						geom_point() + 
						scale_color_manual(values = mrk) + 
						theme(legend.position = "none") + 
						geom_smooth(method="lm",color="black")

pdf(paste(analysis_name,"dat_vol_Da_hatch_non_phylogenetic.pdf",sep=""),width=4,height=4,useDingbats=F)
print(non_phylo_regression)
dev.off()

sink(paste(analysis_name,"non_phylogenetic_regression_volume_time_to_hatching.txt",sep=""))
print(summary(lm(dat_vol_Da_hatch,form=trait2~trait1)))
sink()

# compare egg volume and the interval between mitoses - temperature corrected
egg_vol <- egg_database %>% mutate(trait1 = logvol)
dev_Da_int_btw_preblast_mitoses <- development %>% mutate(trait2 = logDa_int_btw_preblast_mitoses) %>% filter(dipause != "y")
dat_vol_Da_int_btw_preblast_mitoses <- combine_dev_egg_datasets(egg_vol,dev_Da_int_btw_preblast_mitoses)
dat_vol_Da_int_btw_preblast_mitoses_plot <- ggplot(dat_vol_Da_int_btw_preblast_mitoses,aes(x = trait1, y = trait2, color = group)) + 
						geom_point() + 
						scale_color_manual(values = mrk) + 
						theme(legend.position = "none")

pdf(paste(analysis_name,"dat_vol_Da_int_btw_preblast_mitoses.pdf",sep=""),width=4,height=4,useDingbats=F)
print(dat_vol_Da_int_btw_preblast_mitoses_plot)
dev.off()


# compare egg volume and time to midcellularization  - temperature corrected
egg_vol <- egg_database %>% mutate(trait1 = logvol)
dev_Da_time_midcellularization <- development %>% mutate(trait2 = logDa_midcellularization) %>% filter(dipause != "y")
dat_vol_Da_time_midcellularization <- combine_dev_egg_datasets(egg_vol,dev_Da_time_midcellularization)
dat_vol_Da_time_midcellularization_plot <- ggplot(dat_vol_Da_time_midcellularization,aes(x = trait1, y = trait2, color = group)) + 
						geom_point() + 
						scale_color_manual(values = mrk) + 
						theme(legend.position = "none")

pdf(paste(analysis_name,"dat_vol_Da_time_midcellularization.pdf",sep=""),width=4,height=4,useDingbats=F)
print(dat_vol_Da_time_midcellularization_plot)
dev.off()
