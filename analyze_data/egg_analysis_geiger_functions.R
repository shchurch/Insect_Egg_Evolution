### Created by SHC April 2018

### This code provides functions for formatting a data table to analyze with geiger

library(geiger)
library(ape)
library(dplyr)

build_geiger_dataset <- function(data,tree) { 
	data <- data[sample(nrow(data)),] %>% group_by(rank) %>% slice(1L)
	mappers <- data %>% filter(rank %in% tree$tip.label) %>% dplyr::select(rank,trait)
	mps <- na.omit(mappers)
	ms <- mps$trait
	names(ms) <- mps$rank
	return(ms)
}

build_geiger_twotrait_dataset <- function(data,tree) { 
	data <- data[sample(nrow(data)),] %>% group_by(rank) %>% slice(1L)
	mappers <- data %>% filter(rank %in% tree$tip.label) %>% dplyr::select(rank,trait1,trait2)
	mps <- na.omit(mappers)
	ms <- data.frame(trait1 = mps$trait1, trait2 = mps$trait2)
	row.names(ms) <- mps$rank
	return(ms)
}

fit_model <- function(data,tree,model) {
	ms_trim <- drop.tip(tree,setdiff(tree$tip.label,names(data)))
	ms_trim <- multi2di(ms_trim)
	
	fit <- fitContinuous(ms_trim,data,model=model)

	return(fit)
}

