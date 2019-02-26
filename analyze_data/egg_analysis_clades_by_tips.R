### Created by SHC May 2018

### This code executes allometry comparisons across groups of insects determined by number of species

source("analyze_data/egg_analysis_build_dataframe.R")
source("analyze_data/egg_analysis_read_trees.R")
source("analyze_data/egg_analysis_pgls_functions.R")
library(phytools)

### Create a database reduced for the comparison
# Select only the traits of interes
egg_database <- egg_database %>% select(genus,logX1,logX2,group,family,order) %>% 
				# Remove missing data
				na.omit() %>% 
				# Use only tips in the phylogeny
				filter(genus %in% genus_mcc_tree$tip.label)
# Reduce the phylogeny to only tips in the database
tips_not_in_tree <- setdiff(genus_mcc_tree$tip.label,egg_database$genus)
genus_trees <- lapply(genus_trees,drop.tip,tip = tips_not_in_tree)

### Set thresholds for max and min number of species per group
maximum_threshold <- 50
minimum_threshold <- 20

### Determine the groups by finding the minimum number of nodes on the tree that meed the thresholds
# get node numbers
edges <- unique(genus_mcc_tree$edge[,1])
# get tip numbers
tips <- seq(1:(genus_mcc_tree$Nnode + 1))
# get descendants
descendants <- lapply(edges,getDescendants,tree = genus_mcc_tree)
# get number of descendant tips
num_descendants <- sapply(descendants,function(x) {length(x[which(x %in% tips)])})

# find which nodes have fewer than 50 tips
which_sub_50 <- which(num_descendants <= maximum_threshold)
sub_50 <- edges[which_sub_50]
# find the descendants of those nodes
desc_sub_50 <- descendants[which_sub_50] %>% unlist()
# find the smallest set of nodes and tips with no more than 50 descendants
good_nodes <- sub_50[which(!(sub_50 %in% desc_sub_50))]
good_tips <- tips[which(!(tips %in% desc_sub_50))]

# group nodes and tips by sets
sets_node_tips <- c(descendants[match(good_nodes,edges)],good_tips) %>% reshape::melt()
# filter to just tips, order them by the phylogenetic order of a ladderized tree
sets_tips <- sets_node_tips %>% filter(value %in% tips) %>% 
	mutate(tip = genus_mcc_tree$tip.label[value]) %>% 
	arrange(match(tip,genus_mcc_tree$tip.label[genus_mcc_tree$edge[which(genus_mcc_tree$edge[,2] %in% tips),2]]))
# add sets to egg database
egg_database$set <- factor(plyr::mapvalues(egg_database$genus,sets_tips$tip,sets_tips$L1,warn_missing=F),levels=unique(sets_tips$L1))

# determine the groups for subgroup analyses
group_list <- egg_database %>% select(genus,logX1,logX2,set) %>% na.omit() %>% filter(!(is.na(set))) %>% distinct(set,genus) %>% group_by(set) %>% summarise(n = n()) %>% filter(n > minimum_threshold) %>% pull(set)

### Run PGLS comparisons
run_all_allometry_pgls <- function(tree){
	### PGLS length vs width
	length_width <- egg_database %>% select(genus,logX1,logX2,set) %>% rename(rank = genus, trait1 = logX1, trait2= logX2, group = set)
	allometry_pgls_length_width <- run_all_taxa_and_by_group_pgls(length_width,tree,"length_width",group_list)

	return(list(allometry_pgls_length_width))
}

allometry_distribution_raw <- lapply(genus_trees,run_all_allometry_pgls)

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

allometry_group_l_w_table <- lapply(seq(length(group_list)),get_group_allometry_table,pgls=1)

# Print the results
sink("group_clades_by_tips_allometry_tables.txt")
for(i in seq(1:length(group_list))) {
	cat("\n")
	print(as.character(group_list[i]))
	cat("\nlength_width")
	print(allometry_group_l_w_table[[i]])
}
sink()

### Compare the slope of allometric regressions across phylogenetic clades
# This function gets the slopes from a distribution of regression results
get_slope_distribution <- function(pgls,value) {
	slope <- sapply(allometry_distribution_raw,function(x) {x[[pgls]][[2]][[value]]$coefficients[[2]]})
	return(slope)
}

# Get slopes by group
slope_dist_l_w <- lapply(seq(length(group_list)),get_slope_distribution,pgls=1)

# Set up slope plot
sets_groups <- egg_database %>% filter(set %in% group_list) %>% distinct(set,group)
names(slope_dist_l_w) <- group_list
slope_dist_l_w <- reshape::melt(slope_dist_l_w)
slope_dist_l_w$group <- plyr::mapvalues(slope_dist_l_w$L1,sets_groups$set,sets_groups$group)

slope_l_w_plot <- ggplot(slope_dist_l_w,aes(x = rev(L1), y = value, color = group, fill = group)) + 
	geom_abline(intercept=1.0,slope=0,linetype=2) + 
	geom_boxplot(outlier.colour = NULL) + 
	scale_fill_manual(values = mrk) + 
	scale_color_manual(values = mrk) + 
	stat_summary(geom = "crossbar", width=0.65, lwd = 0.7, fatten=0, color="white", fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x))) }) +
	theme(legend.position="none",axis.text.x = element_blank()) + ylab("slope") + xlab("clade")

pdf("slope_l_w_plot_by_clade.pdf",height=3.5,width=7)
print(slope_l_w_plot)
dev.off()

### Create a phylogeny with a single representative per set
single_tips <- egg_database %>% filter(genus %in% genus_mcc_tree$tip.label) %>% filter(set %in% group_list) %>% group_by(set) %>% slice(1L) %>% pull(genus)
reduced_tree <- drop.tip(genus_mcc_tree,setdiff(genus_mcc_tree$tip.label,single_tips))
reduced_tips <- plyr::mapvalues(reduced_tree$tip.label,sets_tips$tip,sets_tips$L1,warn_missing=F)
reduced_tree$tip.label <- paste(reduced_tips,plyr::mapvalues(reduced_tree$tip.label,egg_database$genus,egg_database$order,warn_missing=F),sep="_")

pdf(file="set_phylogeny.pdf")
plot(reduced_tree)
dev.off()

### Create a list of families within each set
set_fam <- egg_database %>% filter(set %in% group_list) %>% distinct(set,family) %>% arrange(set,family)
other_fams <- egg_database %>% filter(!(set %in% group_list)) %>% filter(genus %in% genus_mcc_tree$tip.label) %>% distinct(set,family)
set_fam$unique_to_set <- !(set_fam$family %in% other_fams$family)

write.table(set_fam,file="sets_and_family_list.tsv",sep="\t")

### Create a list of genera within each set
sets_genus_names <- sapply(seq(length(allometry_distribution_raw[[1]][[1]][[2]])),function(x) {names(allometry_distribution_raw[[1]][[1]][[2]][[x]]$fitted)})
names(sets_genus_names) <- group_list

sink(file="sets_genus_names.txt")
print(sets_genus_names)
sink()

save.image(file="egg_analysis_clades_by_tips_workspace.RData")


