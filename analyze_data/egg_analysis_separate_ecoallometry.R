### Created by SHC April 2018

### This code compares the allometric slope across ecological regimes

source("analyze_data/egg_analysis_build_dataframe.R")
source("analyze_data/egg_analysis_read_trees.R")
source("analyze_data/egg_analysis_pgls_functions.R")
library(phytools)
library(corHMM)
library(RColorBrewer)

### Select ecology regime of interest in this script, then source
#source("analyze_data/egg_analysis_parasitoid.R")
source("analyze_data/egg_analysis_aquatic.R")

### select colors
#eco_cols <- c("dark gray","#c55c15")
eco_cols <- c("dark gray","#2056CE")

tree <- genus_mcc_tree

### Select trait of interest
egg_eco_data$trait1 <- egg_eco_data$logX1
egg_eco_data$trait2 <- egg_eco_data$logX2

### Set up the OUwie data frame
reg <- egg_eco_data %>% group_by(rank) %>% slice(1L) %>%
		select(rank,discrete_regime,logX1,logX2,sqcurv,sqasym,family) %>% 
		rename(species = rank,
				discrete = discrete_regime) %>%
		filter(species %in% tree$tip.label)

reg <- as.data.frame(reg)
pp_pruned <- drop.tip(tree,setdiff(tree$tip.label,reg$species))
pp <- rayDISC(pp_pruned,as.data.frame(reg[,c(1,2)]),model="ER",node.states="marginal")

### Plot the ecology ASR on the phylogeny
branch_cols <- c(as.integer(mapvalues(pp$phy$tip.label,reg$species,reg$discrete)),pp$phy$node.label)

pdf(file = paste(analysis_name,".pdf",sep=""),height=12,width=4)
	plot(pp$phy,show.tip.label=F,edge.col=mapvalues(branch_cols[pp$phy$edge[,2]],seq(1:length(eco_cols)),c(eco_cols[-1],"dark gray"))) 
	add.scale.bar()
dev.off()

# get list of edges in the phylogeny
edges <- pp$phy$edge

# create a list of tip states
tip_states <- as.integer(mapvalues(pp$phy$tip.label,reg$species,reg$discrete))
# create a list of node states
node_states <- pp$phy$node.label
# combined list, ordered according to node numbers
states <- c(tip_states,node_states)

# create new state vector
new_states <- rep(0,length(states))

# get starting state
edge_state_1 <- states[edges[,1]]
# get ending state
edge_state_2 <- states[edges[,2]]

# find nodes where start and end state are not the same
edge_diff <- edge_state_2 - edge_state_1
shifts <- edges[which(!(edge_diff == 0)),2]

# get lists of shifts and descendants
shifts_to_state <- edges[which(edge_diff == 1),2]
reversals <- edges[which(edge_diff == -1),2]
reversal_descendants <- lapply(reversals,getDescendants,tree=pp$phy)
num_reversal_descendants <- lapply(reversal_descendants,function(x) length(which(x > length(pp$phy$tip.label))))
reversal_tips <- unlist(lapply(reversal_descendants,function(x){tree$tip.label[x]})) %>% na.omit()

# get descendants of those nodes
shift_descendants <- lapply(shifts,getDescendants,tree=pp$phy)
num_shift_descendants <- lapply(shift_descendants,function(x) length(which(x > length(pp$phy$tip.label))))

# for each of the shifts, set all descendants to a new regime
for (i in (1:length(shift_descendants))) {
	new_states[shift_descendants[[i]]] <- i
	# set the node itself to same regime
	new_states[shifts[i]] <- i
}

# create a new eco regime data frame with separate regimes per shift
new_reg <- data.frame(species = pp$phy$tip.label,shift_regimes = new_states[seq(1,length(pp$phy$tip.label))]) %>% arrange(species)
# new_reg <- left_join(reg,new_reg,by="species") %>% mutate(phy_group = mapvalues(species,egg_database$genus,egg_database$group,warn_missing=F),group_shift_regimes = paste(phy_group,shift_regimes,sep="_"))
new_reg <- left_join(reg,new_reg,by="species") %>% mutate(group_shift_regimes = shift_regimes)

# create a new phylogeny with shift regimes as node labels
new_phy <- pp$phy
new_phy$node.label <- new_states[seq(length(pp$phy$tip.label)+1,length(new_states))]

group_list <- unique(new_reg$group_shift_regimes)

run_all_allometry_pgls <- function(tree){
	### PGLS length vs width
	length_width <- new_reg %>% select(species,logX1,logX2,group_shift_regimes,discrete) %>% rename(rank = species, trait1 = logX1, trait2= logX2, group = group_shift_regimes)
	allometry_pgls_length_width <- run_eco_regime_pgls(length_width,tree,paste(analysis_name,"length_width",sep="_"),group_list)
	
	### PGLS length vs asymmetry, width = independent
	length_asym_width <- new_reg %>% select(species,logX1,logX2,sqasym,group_shift_regimes,discrete) %>% rename(rank = species, trait1 = logX1, trait2= sqasym, indep = logX2, group = group_shift_regimes)
	allometry_pgls_length_asym_width <- run_eco_regime_pgls(length_asym_width,tree,paste(analysis_name,"length_asym_width",sep="_"),group_list)

	### PGLS length vs curvature, width = independent
	length_curv_width <- new_reg %>% select(species,logX1,logX2,sqcurv,group_shift_regimes,discrete) %>% rename(rank = species, trait1 = logX1, trait2= sqcurv, indep = logX2, group = group_shift_regimes)
	allometry_pgls_length_curv_width <- run_eco_regime_pgls(length_curv_width,tree,paste(analysis_name,"length_curv_width",sep="_"),group_list)

	return(list(allometry_pgls_length_width,allometry_pgls_length_asym_width,allometry_pgls_length_curv_width))
}

### Run all PGLS comparisons over the posterior distribution of trees
### Results stored as a list of lists
allometry_distribution_raw <- lapply(genus_trees,run_all_allometry_pgls)
#allometry_distribution_raw <- lapply(list(genus_mcc_tree),run_all_allometry_pgls)

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

### get the genus names used in each pgls comparison
sets_genus_names <- sapply(seq(length(allometry_distribution_raw[[1]][[1]][[2]])),function(x) {names(allometry_distribution_raw[[1]][[1]][[2]][[x]]$fitted)})

sink(file=paste(analysis_name,"sets_genus_names.txt",sep="_"))
print(sets_genus_names)
sink()

groups_l_w <- sapply(seq(length(allometry_distribution_raw[[1]][[1]][[2]])),function(x) {names(allometry_distribution_raw[[1]][[1]][[2]][[x]]$fitted)})
groups_l_asym_w <- sapply(seq(length(allometry_distribution_raw[[1]][[2]][[2]])),function(x) {names(allometry_distribution_raw[[1]][[1]][[2]][[x]]$fitted)})
groups_l_curv_w <- sapply(seq(length(allometry_distribution_raw[[1]][[3]][[2]])),function(x) {names(allometry_distribution_raw[[1]][[3]][[2]][[x]]$fitted)})

allometry_group_l_w_table <- lapply(seq(length(groups_l_w)),get_group_allometry_table,pgls=1)
allometry_group_l_asym_w_table <- lapply(seq(length(groups_l_asym_w)),get_group_allometry_table,pgls=2)
allometry_group_l_curv_w_table <- lapply(seq(length(groups_l_curv_w)),get_group_allometry_table,pgls=3)

sink(paste(analysis_name,"group_clades_by_tips_allometry_tables.txt",sep="_"))
for(i in seq(1:length(allometry_group_l_w_table))) {
	cat("\nlength_width")
	print(allometry_group_l_w_table[[i]])
}
sink()

### Build slope comparisons
get_slope_distribution <- function(pgls,value) {
	slope <- sapply(allometry_distribution_raw,function(x) {x[[pgls]][[2]][[value]]$coefficients[[2]]})
	return(slope)
}

slope_dist_l_w <- lapply(seq(length(groups_l_w)),get_slope_distribution,pgls=1)
groups_ecology <- sapply(seq(length(groups_l_w)),function(x) {new_reg %>% filter(species == groups_l_w[[x]][1]) %>% pull(discrete)})
groups_regimes <- sapply(seq(length(groups_l_w)),function(x) {new_reg %>% filter(species == groups_l_w[[x]][1]) %>% pull(group_shift_regimes)})

names(slope_dist_l_w) <- groups_regimes
slope_dist_l_w <- reshape::melt(slope_dist_l_w)
slope_dist_l_w$ecology <- mapvalues(slope_dist_l_w$L1,groups_regimes,groups_ecology)

slope_l_w_plot <- ggplot(slope_dist_l_w,aes(x = L1, y = value, color = ecology, fill = ecology)) + 
	geom_abline(intercept=1.0,slope=0,linetype=2) + 
	geom_boxplot(outlier.colour = NULL) + 
	scale_fill_manual(values = eco_cols) + 
	scale_color_manual(values = eco_cols) + 
	stat_summary(geom = "crossbar", width=0.65, lwd = 0.7, fatten=0, color="white", fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x))) }) +
	theme(legend.position="none",axis.text.x = element_blank()) + ylab("slope") + xlab("clade")

pdf(paste(analysis_name,"l_w_separate_ecoallometry_slope_plot.pdf",sep="_"),height=4,width=7)
print(slope_l_w_plot)
dev.off()

save.image(file=paste(analysis_name,"separate_eco_allometry_workspace.RData",sep="_"))

