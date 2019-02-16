### Created by SHC April 2018

### This code reads in the posterior and MCC trees for phylogenetic analyses

library(ape)
library(phangorn)
library(phytools)
library(dplyr)

## Read in distribution of ultrametric trees
## Select the backbone - Misof or Rainford
genus_trees <- read.nexus("phylogeny/final_trees/misof_posterior_ultrametric.nxs")
rainford_genus_trees <- read.nexus("phylogeny/final_trees/rainford_posterior_ultrametric.nxs")

### Read in the mcc tree
genus_mcc_tree <- read.nexus("phylogeny/final_trees/misof_mcc_ultrametric.nxs")
rainford_genus_mcc_tree <- read.nexus("phylogeny/final_trees/rainford_mcc_ultrametric.nxs")

## Read the family level tree from the Rainford publication
fam_tree <- read.nexus("phylogeny/final_trees/rainford_family_ultrametric.nxs")

## Plot morphospace covered by tree
##fam_morphospace_plot <- ggplot(egg_database,aes(x = logar, y = logvol)) + geom_point(color = "dark grey",alpha=0.5,stroke=0,size=1) + geom_point(data = egg_database %>% filter(family %in% fam_tree$tip.label),color = "blue",alpha=0.5,stroke=0,size=1)
##pdf(file="fam_morphospace_plot.pdf",width=4,height=4)
##	print(fam_morphospace_plot)
##dev.off()

##genus_morphospace_plot <- ggplot(egg_database,aes(x = logar, y = logvol)) + geom_point(color = "dark grey",alpha=0.5,stroke=0,size=1) + geom_point(data = egg_database %>% filter(genus %in% genus_mcc_tree$tip.label),color = "red",alpha=0.5,stroke=0,size=1)
##pdf(file="genus_morphospace_plot.pdf",width=4,height=4)
##	print(genus_morphospace_plot)
##dev.off()

## Break the tips into discrete monophyletic clades
## determined by number of tips

# get node numbers
edges <- unique(genus_mcc_tree$edge[,1])
# get tip numbers
tips <- seq(1:(genus_mcc_tree$Nnode + 1))
# get descendants
descendants <- lapply(edges,getDescendants,tree = genus_mcc_tree)
# get number of descendant tips
num_descendants <- sapply(descendants,function(x) {length(x[which(x %in% tips)])})
# find which nodes have fewer than 50 tips
which_sub_50 <- which(num_descendants <= 50)
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

### These commands were used once to clean up and format the trees

### MISOF BACKBONE
### read the distribution of genus level trees from our phylogenetic analysis
#	genus_trees <- read.nexus("phylogeny/final_trees/misof_100_posterior.trees")
#	### ladderize them for easy visualization
#	genus_trees <- lapply(genus_trees,ladderize)
#	### sub out the spaces in tip labels
#	for(i in 1:length(genus_trees)){
#		genus_trees[[i]]$tip.label <- gsub("_[0-9]+","",genus_trees[[i]]$tip.label)
#	
#		### force family tree to be ultrametric
#		genus_trees[[i]]<-nnls.tree(cophenetic(genus_trees[[i]]),genus_trees[[i]],rooted=TRUE)
#	}
#	### write ultrametric trees
#	write.nexus(file="phylogeny/final_trees/misof_posterior_ultrametric.nxs",genus_trees)
#	
#	### RAINFORD BACKBONE
#	### read the distribution of genus level trees from our phylogenetic analysis
#	genus_trees <- read.nexus("phylogeny/final_trees/rainford_100_posterior.trees")
#	### ladderize them for easy visualization
#	genus_trees <- lapply(genus_trees,ladderize)
#	### sub out the spaces in tip labels
#	for(i in 1:length(genus_trees)){
#		genus_trees[[i]]$tip.label <- gsub("_[0-9]+","",genus_trees[[i]]$tip.label)
#	
#		### force family tree to be ultrametric
#		genus_trees[[i]]<-nnls.tree(cophenetic(genus_trees[[i]]),genus_trees[[i]],rooted=TRUE)
#	}
#	### write ultrametric trees
#	write.nexus(file="phylogeny/final_trees/rainford_posterior_ultrametric.nxs",genus_trees)

#	### MISOF MCC
#	### read the distribution of genus level trees from our phylogenetic analysis
#	genus_mcc_tree <- read.nexus("phylogeny/final_trees/misof_mcc.tree")
#	### ladderize them for easy visualization
#	genus_mcc_tree <- ladderize(genus_mcc_tree)
#	### sub out the spaces in tip labels
#	genus_mcc_tree$tip.label <- gsub("_[0-9]+","",genus_mcc_tree$tip.label)
#	### force family tree to be ultrametric
#	genus_mcc_tree<-nnls.tree(cophenetic(genus_mcc_tree),genus_mcc_tree,rooted=TRUE)
#	### write ultrametric trees
#	write.nexus(file="phylogeny/final_trees/misof_mcc_ultrametric.nxs",genus_mcc_tree)
#
#	### RAINFORD MCC
#	### read the distribution of genus level trees from our phylogenetic analysis
#	genus_mcc_tree <- read.nexus("phylogeny/final_trees/rainford_mcc.tree")
#	### ladderize them for easy visualization
#	genus_mcc_tree <- ladderize(genus_mcc_tree)
#	### sub out the spaces in tip labels
#	genus_mcc_tree$tip.label <- gsub("_[0-9]+","",genus_mcc_tree$tip.label)
#	### force family tree to be ultrametric
#	genus_mcc_tree<-nnls.tree(cophenetic(genus_mcc_tree),genus_mcc_tree,rooted=TRUE)
#	### write ultrametric trees
#	write.nexus(file="phylogeny/final_trees/rainford_mcc_ultrametric.nxs",genus_mcc_tree)
#	
#	### RAINFORD FAMILY TREE
#	fam_trees <- read.nexus("phylogeny/final_trees/rainford_family_tree.nxs")
#	fam_tree <- fam_trees[[2]]
#	### Get tip reference table from Rainford publication
#	fam_tree_key <- read.delim("phylogeny/final_trees/rainford_family_key.txt",header=F,stringsAsFactors=F)
#	### Change tip labels to families / orders
#	fam_tree$tip.label <- plyr::mapvalues(fam_tree$tip.label,fam_tree_key$V3,fam_tree_key$V2)
#	### Force family tree to be ultrametric
#	fam_tree<-nnls.tree(cophenetic(fam_tree),fam_tree,rooted=TRUE)
#	### write ultrametric trees
#	write.nexus(file="phylogeny/final_trees/rainford_family_ultrametric.nxs",fam_tree)


