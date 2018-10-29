### This file was created by SCH April 2018

### This code compares OU models of egg size / shape with ecological regimes

library(corHMM)
library(OUwie)

### Select ecology regime of interest in this script, then source
#source("analyze_data/egg_analysis_parasitoid.R")
source("analyze_data/egg_analysis_aquatic.R")
#source("analyze_data/egg_analysis_ovary.R")
#source("analyze_data/egg_analysis_wingless_phasmatodea.R")
#source("analyze_data/egg_analysis_migratory_lepidoptera.R")

analysis_name <- paste("ouwie",analysis_name,sep="_")
### select colors
#eco_cols <- c("dark gray","#c55c15")
eco_cols <- c("dark gray","#2056CE")
#eco_cols <- c("dark gray","#2056CE","#c55c15")

### Select trait of interest
#egg_eco_data$trait1 <- egg_eco_data$logvol
#analysis_name <- paste("logvol",analysis_name,sep="_")

#egg_eco_data$trait1 <- egg_eco_data$logar
#analysis_name <- paste("logar",analysis_name,sep="_")

#egg_eco_data$trait1 <- egg_eco_data$sqcurv
#analysis_name <- paste("sqcurv",analysis_name,sep="_")

egg_eco_data$trait1 <- egg_eco_data$sqasym
analysis_name <- paste("sqasym",analysis_name,sep="_")

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
branch_cols <- c(as.integer(mapvalues(pp$phy$tip.label,reg$species,reg$discrete)),pp$phy$node.label)

pdf(file = paste(analysis_name,".pdf",sep=""),height=12,width=4)
	plot(pp$phy,show.tip.label=F,edge.col=mapvalues(branch_cols[pp$phy$edge[,2]],seq(1:length(eco_cols)),c(eco_cols[-1],"dark gray"))) 
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
