### Created by SHC  May 2018

### These functions are used to iteratively annotate the egg database by taxonomic groups

### This function accepts binary states (coded as TRUE or FALSE)
update_regimes_TF <- function(eco_rank,egg_names){
	# get taxa that are in each state
	eco_true <- eco_data %>% filter(rank == eco_rank) %>% filter(ecology == "TRUE") %>% pull(name)
	eco_false <- eco_data %>% filter(rank == eco_rank) %>% filter(ecology == "FALSE") %>% pull(name)
	eco_some <- eco_data %>% filter(rank == eco_rank) %>% filter(ecology == "SOME") %>% pull(name)

	# update database
	if(length(which(egg_names %in% eco_true)) > 0){
		egg_eco_data[which(egg_names %in% eco_true),]$eco_regime <- "derived"
	}
	if(length(which(egg_names %in% eco_false)) > 0){
		egg_eco_data[which(egg_names %in% eco_false),]$eco_regime <- "ancestral"
	}
	#### CURRENTLY SETS SOME TO TRUE
	if(length(which(egg_names %in% eco_some)) > 0){
		egg_eco_data[which(egg_names %in% eco_some),]$eco_regime <- "derived"
		#egg_eco_data[which(egg_names %in% eco_some),]$eco_regime <- "ancestral"
	}

	return(egg_eco_data)
}

### This function accepts multistates, encoded in eco_regimes 
update_regimes_multistate <- function(eco_rank,egg_names,eco_regimes){
	# get taxa that are in each state
	eco_1 <- eco_data %>% filter(rank == eco_rank) %>% filter(ecology == eco_regimes[[1]]) %>% pull(name)
	eco_2 <- eco_data %>% filter(rank == eco_rank) %>% filter(ecology == eco_regimes[[2]]) %>% pull(name)
	eco_3 <- eco_data %>% filter(rank == eco_rank) %>% filter(ecology == eco_regimes[[3]]) %>% pull(name)
	eco_4 <- eco_data %>% filter(rank == eco_rank) %>% filter(ecology == eco_regimes[[4]]) %>% pull(name)

	# update database 
	if(length(which(egg_names %in% eco_1)) > 0){
		egg_eco_data[which(egg_names %in% eco_1),]$eco_regime <- eco_regimes[[1]]
	}
	if(length(which(egg_names %in% eco_2)) > 0){
		egg_eco_data[which(egg_names %in% eco_2),]$eco_regime <- eco_regimes[[2]]
	}
	if(length(which(egg_names %in% eco_3)) > 0){
		egg_eco_data[which(egg_names %in% eco_3),]$eco_regime <- eco_regimes[[3]]
	}
	if(length(which(egg_names %in% eco_4)) > 0){
		egg_eco_data[which(egg_names %in% eco_4),]$eco_regime <- eco_regimes[[4]]
	}

	return(egg_eco_data)
}
