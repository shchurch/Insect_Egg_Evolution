### This file was created by SCH April 2018

### This code compares OU models of egg size / shape with ecological regimes
source("analyze_data/egg_analysis_build_dataframe.R")
source("analyze_data/egg_analysis_read_trees.R")
library(corHMM)
library(OUwie)

args = commandArgs(trailingOnly=TRUE)
### This analysis takes 5 arguments
# 1 = 'relaxed' or 'strict' ecological classification method
# 2 = using the Misf 'mcc' tree, 'Rainford_mcc', or an integer N for the Misof posterior distribution tree
# 3 = using 'observed' or 'simulated' ecology states 
# 4 = testing 'parasitoid', 'internal', 'aquatic', 'combined_aquatic', 'wingless_Phasmatodea',  'migratory_Lepidoptera', or 'three_state'
# 5 = comparing the trait 'vol', 'ar', 'asym', or 'curv', all transformed 

analysis_name <- "ouwie"

### Select classificaiton method
class_flag <- args[1]
analysis_name <- paste(analysis_name,args[1],sep="_")

### Select tree 
# Misof maximum clade credibility tree
if(args[2] == "mcc") {
	tree <- genus_mcc_tree
	analysis_name <- paste(analysis_name,"mcc",sep="_")
# Rainfored mcc tree
} else if(args[2] == "mcc_Rainford") {
	tree <- rainford_genus_mcc_tree
	analysis_name <- paste(analysis_name,"mcc_rainford",sep="_")
# Posterior distribution
} else {
	# set random seed for simulation to the integer of the array
	set.seed(as.numeric(args[2]))
	array_num <- as.numeric(args[2])
	tree <- genus_trees[[array_num]]
	analysis_name <- paste(analysis_name,"posterior",sep="_")
	analysis_name <- paste(analysis_name,array_num,sep="_")
}

### Select ecology regime of interest in this script
ecology_arg <- args[4]
if(args[3] == "observed") {
	if(ecology_arg %in%  c("internal", "parasitoid")) {
		source("analyze_data/egg_analysis_parasitoid.R")

	} else if(ecology_arg %in% c("aquatic", "combined_aquatic")) {
		source("analyze_data/egg_analysis_aquatic.R")

	} else if(ecology_arg == "wingless_Phasmatodea") {
		source("analyze_data/egg_analysis_wingless_phasmatodea.R")

	} else if(ecology_arg == "migratory_Lepidoptera") {
		source("analyze_data/egg_analysis_migratory_lepidoptera.R")

	} else if(ecology_arg == "three_state") {
		ecology_arg <- "internal"
		source("analyze_data/egg_analysis_parasitoid.R")
		eco_data_internal <- egg_eco_data

		ecology_arg <- "aquatic"
		source("analyze_data/egg_analysis_aquatic.R")
		eco_data_aquatic <- egg_eco_data

		egg_eco_data <- eco_data_internal %>% rename(internal = discrete_regime) %>% 
			left_join((eco_data_aquatic %>% 
				select(rank,discrete_regime) %>% 
				rename(aquatic = discrete_regime)),by="rank") %>% 
			mutate(discrete_regime = ifelse(internal == 2,internal,ifelse(aquatic == 2, 3,1)))
	
		analysis_name <- paste(analysis_name,"three_state",sep="_")

		eco_cols <- c("dark gray","#c55c15","#2056CE")
	} else {
		warning("incorrect ecology argument")
		quit(status=1)
	}
} else if(args[3] == "simulated") {
	source("analyze_data/egg_analysis_simulate_ecology.R")
	analysis_name <- paste(analysis_name,"simulated",sep="_")

} else {
	warning("incorrect observed / simulated argument")
	quit(status=1)
}

### Select trait of interest
# volume
if(args[5] == "vol") {
	egg_eco_data$trait1 <- egg_eco_data$logvol
	analysis_name <- paste(analysis_name,"logvol",sep="_")
# aspect ratio
} else if(args[5] == "ar") {
	egg_eco_data$trait1 <- egg_eco_data$logar
	analysis_name <- paste(analysis_name,"logar",sep="_")
# curvature
} else if(args[5] == "curv") {
	egg_eco_data$trait1 <- egg_eco_data$sqcurv
	analysis_name <- paste(analysis_name,"sqcurv",sep="_")
# asymmetry
} else if(args[5] == "asym") {
	egg_eco_data$trait1 <- egg_eco_data$sqasym
	analysis_name <- paste(analysis_name,"sqasym",sep="_")

} else {
	warning("incorrect trait argument")
	quit(status=1)
}

### Set up the OUwie data frame
reg <- egg_eco_data %>% 
		# give variables the names expected by OUwie
		select(rank,discrete_regime,trait1) %>% 
		rename(species = rank,
				discrete = discrete_regime,
				continuous = trait1) %>%
		# remove NA
		na.omit() %>% 
		# (optional) maximize number in derived state
		arrange(species,desc(discrete)) %>% 
		# randomly choose one representative per genus
		group_by(species) %>% slice(1L) %>%
		filter(species %in% tree$tip.label) %>% as.data.frame()

### Prune the full tree to only those with ecological + trait data
pp_pruned <- drop.tip(tree,setdiff(tree$tip.label,reg$species))

### corHMM, using the rayDISC function, reconstructs ancestral states
pp <- rayDISC(pp_pruned,as.data.frame(reg[,c(1,2)]),model="ER",node.states="marginal")
 	
### Plot the ecology ASR on the phylogeny
branch_cols <- c(as.integer(plyr::mapvalues(pp$phy$tip.label,reg$species,reg$discrete)),pp$phy$node.label)

pdf(file = paste(analysis_name,".pdf",sep=""),height=12,width=4)
	plot(pp$phy,show.tip.label=F,edge.col=plyr::mapvalues(branch_cols[pp$phy$edge[,2]],seq(1:length(eco_cols)),eco_cols)) 
	add.scale.bar()
dev.off()

### Run with 3 models
ou1_res <- OUwie(pp$phy,reg,"OU1",simmap.tree=FALSE,diagn=FALSE)
oum_res <- OUwie(pp$phy,reg,"OUM",simmap.tree=FALSE,diagn=FALSE)
bm1_res <- OUwie(pp$phy,reg,"BM1",simmap.tree=FALSE,diagn=FALSE)

### Print the results
results <- list(bm1_res,ou1_res,oum_res)
names(results) <- c("BM1","OU1","OUM")

results_table <- data.frame(row.names = c("Brownian Motion", "OU - 1 Optimum", "OU - Multiple Optima"), "AICc" = c(results[[1]]$AICc,results[[2]]$AICc,results[[3]]$AICc))

sink(paste(analysis_name,"_results_table.txt",sep=""))
print(results_table)
sink()

save.image(file = paste(analysis_name,".Rdata",sep=""))
