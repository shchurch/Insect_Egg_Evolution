library(dplyr)

results_dir <- "results_Feb2019/ouwie_results_Feb2019/"

get_aicc <- function(filename) {
	f <- paste(results_dir,filename,".txt",collapse="",sep="")
	lines <- readLines(f)
	aicc <- unlist(lapply(lines[2:4],stringr::str_extract,"-*[0-9]+\\.*[0-9]*$"))
	aicc_df <- data.frame(stringsAsFactors=F,analysis = c("BM1","OU1","OUM"),AICc = as.numeric(aicc))
	return(aicc_df)
}

get_delta_aicc <- function(a) {
	d_a <- a %>% mutate(deltaAICc = AICc - min(AICc))
	return(d_a)
}

get_count_aicc <- function(a) { 
	oum_bm1 <- a$AICc[1] - a$AICc[3]
	oum_ou1 <- a$AICc[2] - a$AICc[3]
	ou1_bm1 <- a$AICc[1] - a$AICc[2]
	res <- FALSE
	if(oum_bm1 > 2 && oum_ou1 > 2){
		res <- "OUM"
	} else if (ou1_bm1 > 2) {
		res <- "OU1"
	} else {
		res <- "BM1"
	}
	return(res)
}

get_filenames <- function(class,tree,eco,data,trait) {
	filenames <- NULL
	if(data == "simulated") {
		filenames <- apply(expand.grid("ouwie",class,tree,as.character(seq(1:100)),eco,data,trait,"results_table"), 1, paste, collapse="_")
	} else {
		filenames <- apply(expand.grid("ouwie",class,tree,as.character(seq(1:100)),eco,trait,"results_table"), 1, paste, collapse="_")
	}
	return(as.list(filenames))
}

get_bm1_delta <- function(a) {
	bm1 <- a$AICc[1] - a$AICc[3]
	return(bm1)
}

get_ou1_delta <- function(a) {
	ou1 <- a$AICc[2] - a$AICc[3]
	return(ou1)
}

get_posterior_results <- function(trait,class,data,eco,tree) {
	analysis_names_df <- expand.grid(class,tree,eco,data,trait,stringsAsFactors=F)
	filenames <- apply(analysis_names_df,1,function(x)get_filenames(x[1],x[2],x[3],x[4],x[5]))

	aicc <- lapply(filenames,lapply,get_aicc)
	delta_aicc <- lapply(aicc,lapply,get_delta_aicc)
	count_delta_aicc <- lapply(aicc,function(x) {table(unlist(lapply(x,get_count_aicc)))})

	avg_bm1_delta <- lapply(aicc,function(x) {mean(unlist(lapply(x,get_bm1_delta)))})
	avg_ou1_delta <- lapply(aicc,function(x) {mean(unlist(lapply(x,get_ou1_delta)))})

	count_delta_aicc <- lapply(aicc,function(x) {table(unlist(lapply(x,get_count_aicc)))['TRUE']})

	results <- data.frame(analysis = apply(analysis_names_df,1,paste,collapse="_"),
			OUM_best_count = unlist(count_delta_aicc),
			avg_BM1_OUM_delta = unlist(avg_bm1_delta),
			avg_OU1_OUM_delta = unlist(avg_ou1_delta)) %>% 
			mutate(OUM_best_count = ifelse(!is.na(OUM_best_count),OUM_best_count,0),
				OUM_best_freq = (OUM_best_count / 100))
	return(results)
}

get_simulated_results <- function(trait,class,data,eco,tree) {
	analysis_names_df <- expand.grid(class,tree,eco,data,trait,stringsAsFactors=F)
	filenames <- apply(analysis_names_df,1,function(x)get_filenames(x[1],x[2],x[3],x[4],x[5]))

	aicc <- lapply(filenames,lapply,get_aicc)
	delta_aicc <- lapply(aicc,lapply,get_delta_aicc)
	count_delta_aicc <- lapply(aicc,function(x) {table(unlist(lapply(x,get_count_aicc)))})

	avg_bm1_delta <- lapply(aicc,function(x) {mean(unlist(lapply(x,get_bm1_delta)))})
	avg_ou1_delta <- lapply(aicc,function(x) {mean(unlist(lapply(x,get_ou1_delta)))})

	min_bm1_delta <- lapply(aicc,function(x) {min(unlist(lapply(x,get_bm1_delta)))})
	min_ou1_delta <- lapply(aicc,function(x) {min(unlist(lapply(x,get_ou1_delta)))})

	bm1_delta <- lapply(aicc,lapply,get_bm1_delta)
	ou1_delta <- lapply(aicc,lapply,get_ou1_delta)

	count_delta_aicc <- lapply(aicc,function(x) {table(unlist(lapply(x,get_count_aicc)))['TRUE']})

	joint_pval <- data.frame(ou1_delta = unlist(ou1_delta[[4]]),
						bm1_delta = unlist(bm1_delta[[4]]),
						min_bm1 = min_bm1_delta[[3]], 
						min_ou1 = min_ou1_delta[[3]]) %>% 
				mutate(ou1_true = ifelse(ou1_delta > min_ou1,'TRUE','FALSE'),
					bm1_true = ifelse(bm1_delta > min_bm1,'TRUE','FALSE'),
					both_true = ifelse(bm1_true == 'TRUE' & ou1_true == 'TRUE','TRUE','FALSE')) %>% select(both_true) %>% length()


	results <- data.frame(analysis = apply(analysis_names_df,1,paste,collapse="_"),
			OUM_best_count = unlist(count_delta_aicc)) %>% 
			mutate(OUM_best_count = ifelse(!is.na(OUM_best_count),OUM_best_count,0),
				OUM_best_freq = (OUM_best_count / 100),
				BM1_boot_pval = c(NA,
					length(which(unlist(bm1_delta[[2]]) > min_bm1_delta[[1]])) / 100,
					NA,
					length(which(unlist(bm1_delta[[4]]) > min_bm1_delta[[3]])) / 100),
				OU1_boot_pval = c(NA,
					length(which(unlist(ou1_delta[[2]]) > min_ou1_delta[[1]])) / 100,
					NA,
					length(which(unlist(ou1_delta[[4]]) > min_ou1_delta[[3]])) / 100),
				joint_pval = c(NA,
					(data.frame(ou1_delta = unlist(ou1_delta[[2]]),
						bm1_delta = unlist(bm1_delta[[2]]),
						min_bm1 = min_bm1_delta[[1]], 
						min_ou1 = min_ou1_delta[[1]]) %>% 
					mutate(ou1_true = ifelse(ou1_delta > min_ou1,'TRUE','FALSE'),
						bm1_true = ifelse(bm1_delta > min_bm1,'TRUE','FALSE'),
						both_true = ifelse(bm1_true == 'TRUE' & ou1_true == 'TRUE','TRUE','FALSE')) %>% 
					select(both_true) %>% length()) / 100,
					NA,
					(data.frame(ou1_delta = unlist(ou1_delta[[4]]),
							bm1_delta = unlist(bm1_delta[[4]]),
							min_bm1 = min_bm1_delta[[3]], 
							min_ou1 = min_ou1_delta[[3]]) %>% 
					mutate(ou1_true = ifelse(ou1_delta > min_ou1,'TRUE','FALSE'),
						bm1_true = ifelse(bm1_delta > min_bm1,'TRUE','FALSE'),
						both_true = ifelse(bm1_true == 'TRUE' & ou1_true == 'TRUE','TRUE','FALSE')) %>% 
					select(both_true) %>% length()) / 100))
	return(results)
}

get_mcc_results <- function(trait,class,data,eco,tree) {
	filenames <- apply(expand.grid("ouwie",class,tree,eco,trait,"results_table"), 1, paste, collapse="_")
	aicc <- lapply(filenames,get_aicc)
	delta_aicc <- lapply(aicc,get_delta_aicc)
	count_delta_aicc <- sapply(aicc,get_count_aicc)
		
	results <- data.frame(analysis = filenames,
			BM1_OUM_delta = sapply(aicc,get_bm1_delta),
			OU1_OUM_delta = sapply(aicc,get_ou1_delta),
			OUM_best = sapply(aicc,get_count_aicc))
	return(results)
}

## MAIN
main_results <-  get_posterior_results(trait = c("logvol","logar","sqasym","sqcurv"),
	class = c("relaxed"),
	data = c("observed"),
	eco = c("internal","in_water","migratory_Lepidoptera","wingless_Phasmatodea"),
	tree = c("posterior"))

## simulated internal 
simulated_internal_results <- get_simulated_results(trait = c("logvol","sqasym"),
	class = c("relaxed"),
	data = c("observed","simulated"),
	eco = c("internal"),
	tree = c("posterior"))

## simulated in water 
simulated_in_water_results <- get_simulated_results(trait = c("logvol","logar"),
	class = c("relaxed"),
	data = c("observed","simulated"),
	eco = c("in_water"),
	tree = c("posterior"))

## mcc_rainford
rainford_mcc_results <-  get_mcc_results(trait = c("logvol","logar","sqasym","sqcurv"),
	class = c("relaxed"),
	data = c("observed"),
	eco = c("internal","in_water"),
	tree = c("mcc_rainford"))

## mcc_misof
misof_mcc_results <-  get_mcc_results(trait = c("logvol","logar","sqasym","sqcurv"),
	class = c("relaxed"),
	data = c("observed"),
	eco = c("internal","in_water"),
	tree = c("mcc"))

## classification method
strict_results <-  get_mcc_results(trait = c("logvol","logar","sqasym","sqcurv"),
	class = c("strict"),
	data = c("observed"),
	eco = c("internal","in_water"),
	tree = c("mcc"))

## ecological definition results
eco_definition_results <-  get_mcc_results(trait = c("logvol","logar","sqasym","sqcurv"),
	class = c("relaxed"),
	data = c("observed"),
	eco = c("parasitoid","combined_aquatic"),
	tree = c("mcc"))

