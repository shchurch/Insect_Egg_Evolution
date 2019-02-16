### This file was created by SHC April 2018

### This code compares OU models of egg size / shape with ecological regimes using OUCH
source("analyze_data/egg_analysis_build_dataframe.R")
source("analyze_data/egg_analysis_read_trees.R")
library(corHMM)
library(ouch)

args = commandArgs(trailingOnly=TRUE)
### This analysis takes 4 arguments
# 1 = 'relaxed' or 'strict' ecological classification method
# 2 = using the Misf 'mcc' tree, 'Rainford_mcc', or an integer N for the Misof posterior distribution tree
# 3 = using 'observed' or 'simulated' ecology states 
# 4 = testing 'parasitoid', 'internal', 'aquatic', 'combined_aquatic', 'wingless_Phasmatodea',  'migratory_Lepidoptera', or 'three_state'
# args[1] <- "relaxed"
# args[2] <- "mcc"
# args[3] <- "observed"
# args[4] <- "internal"

analysis_name <- "ouch"

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

# Select the traits for analysis
egg_eco_data$trait1 <- egg_eco_data$logX1 #length
egg_eco_data$trait2 <- egg_eco_data$logX2 #width
egg_eco_data$trait3 <- egg_eco_data$sqasym #asymmetry
egg_eco_data$trait4 <- egg_eco_data$sqcurv #curvature
analysis_name <- paste(analysis_name,"multivariate",sep="_")

### Set up the OUwie data frame
reg <- egg_eco_data %>% 
		# give variables the names expected by OUwie
		select(rank,discrete_regime,trait1,trait2,trait3,trait4) %>% 
		rename(species = rank,
				discrete = discrete_regime) %>%
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

### Build the dataframe for ouch
ou_tree <- ape2ouch(pp$phy)
ou_tree_df <- as(ou_tree,"data.frame")
ou_reg <- reg %>% rename(labels = species)
ou_data <- merge(ou_tree_df,ou_reg,by="labels",all=T)
row.names(ou_data) <- ou_data$nodes
new_ou_tree <- ouchtree(nodes=ou_data$nodes,
				ancestors = ou_data$ancestors,
				times = ou_data$times,
				labels = ou_data$labels)

### Brownian Motion
brown_fit <- brown(data=ou_data[c("trait1","trait2","trait3","trait4")],tree=new_ou_tree)

### OU - single regime
ou_data$single <- as.factor("single")
ou1_fit <- hansen(data=ou_data[c("trait1","trait2","trait3","trait4")],tree=new_ou_tree, regimes=ou_data["single"], sqrt.alpha=c(1,0,0,0,1,0,0,1,0,1), sigma=c(1,0,0,0,1,0,0,1,0,1), fit=TRUE)

### OU - multiple regimes
ou_data$discrete <- plyr::mapvalues(ou_data$labels,ou_reg$labels,ou_reg$discrete)
oum_fit <- hansen(data=ou_data[c("trait1","trait2","trait3","trait4")],tree=new_ou_tree, regimes=ou_data["discrete"], sqrt.alpha=c(1,0,0,0,1,0,0,1,0,1), sigma=c(1,0,0,0,1,0,0,1,0,1), fit=TRUE)

### Print results
sink(paste(analysis_name,"_results_table.txt",sep=""))
print(summary(brown_fit))
print(summary(ou1_fit))
print(summary(oum_fit))
sink()

save.image(file = paste(analysis_name,".Rdata",sep=""))