### Created by SHC April 2018

### This code compares the fit of evolutionary models for egg size and shape parameters

source("analyze_data/egg_analysis_build_dataframe.R")
source("analyze_data/egg_analysis_read_trees.R")
source("analyze_data/egg_analysis_geiger_functions.R")
library(geiger)
library(arbutus)
library(phytools)

fit_models <- function(data,tree) {
	pruned <- drop.tip(tree,setdiff(tree$tip.label,names(data)))
	
	# fit a Brownian Motion model
	B <- fitContinuous(pruned,data,model="BM")
	# fit an Ornstein-Uhlenbeck model
	O <- fitContinuous(pruned,data,model="OU")
	# fit an Early Burst model
	E <- fitContinuous(pruned,data,model="EB")
	# fit a non phylogenetic model
	W <- fitContinuous(pruned,data,model="white")
	
	fits <- list(B,O,E,W)
	return(fits)
}

### Use the MCC summary tree
genus_tree <- genus_mcc_tree
egg_fit_data <- egg_database %>% mutate(rank = genus)

### Build body size dataset
# Set factor for resampling the body size dataframe, here no downsampling used
downsample_factor <- 1.0
fam_count_threshold <- 0
source("analyze_data/egg_analysis_body_size.R")

### Build the geiger datasets for each egg size and shape parameter
vol_fit_data  <- egg_fit_data %>% mutate(trait = logvol) %>% build_geiger_dataset(tree = genus_tree)
ar_fit_data   <- egg_fit_data %>% mutate(trait = logar) %>% build_geiger_dataset(tree = genus_tree)
asym_fit_data <- egg_fit_data %>% mutate(trait = sqasym) %>% build_geiger_dataset(tree = genus_tree)
curv_fit_data <- egg_fit_data %>% mutate(trait = sqcurv) %>% build_geiger_dataset(tree = genus_tree)
X1_fit_data   <- egg_fit_data %>% mutate(trait = logX1) %>% build_geiger_dataset(tree = genus_tree)
X2_fit_data   <- egg_fit_data %>% mutate(trait = logX2) %>% build_geiger_dataset(tree = genus_tree)

### Build the geiger datasets for body size parameters
body_fit_data <- egg_database_family_body %>% mutate(rank = family, trait = logbodyvol) %>% build_geiger_dataset(tree = fam_tree)
vol_fam_fit_data <- egg_database_family_body %>% mutate(rank = family, trait = logvol) %>% build_geiger_dataset(tree = fam_tree)

### Fit evolutionary models for each egg size, egg shape, and body size parameter
vol_fit  <- vol_fit_data[order(match(names(vol_fit_data),genus_tree$tip.label))] %>% fit_models(tree = genus_tree)
ar_fit   <- ar_fit_data[order(match(names(ar_fit_data),genus_tree$tip.label))] %>% fit_models(tree = genus_tree)
asym_fit <- asym_fit_data[order(match(names(asym_fit_data),genus_tree$tip.label))] %>% fit_models(tree = genus_tree)
curv_fit <- curv_fit_data[order(match(names(curv_fit_data),genus_tree$tip.label))] %>% fit_models(tree = genus_tree)
X1_fit   <- X1_fit_data[order(match(names(X1_fit_data),genus_tree$tip.label))] %>% fit_models(tree = genus_tree)
X2_fit   <- X2_fit_data[order(match(names(X2_fit_data),genus_tree$tip.label))] %>% fit_models(tree = genus_tree)
body_fit <- body_fit_data[order(match(names(body_fit_data),fam_tree$tip.label))] %>% fit_models(tree = fam_tree)
vol_fam_fit <- vol_fam_fit_data[order(match(names(vol_fam_fit_data),fam_tree$tip.label))] %>% fit_models(tree = fam_tree)

### Calculate the significance of phylogenetic signal for each parameter
vol_lambda  <- phylosig(genus_tree,vol_fit_data,test=T,method="lambda")
ar_lambda   <- phylosig(genus_tree,ar_fit_data,test=T,method="lambda")
asym_lambda <- phylosig(genus_tree,asym_fit_data,test=T,method="lambda")
curv_lambda <- phylosig(genus_tree,curv_fit_data,test=T,method="lambda")
X1_lambda   <- phylosig(genus_tree,X1_fit_data,test=T,method="lambda")
X2_lambda   <- phylosig(genus_tree,X2_fit_data,test=T,method="lambda")
body_lambda <- phylosig(fam_tree,body_fit_data,test=T,method="lambda")
vol_fam_lambda <- phylosig(fam_tree,vol_fam_fit_data,test=T,method="lambda")

### Create a table summarizing the results
fit_table <- data.frame(row.names = c("Volume","Aspect Ratio","Asymmetry","Curvature","Length","Width","Cubic body length - family","Egg volume - family"),
	BM = c(vol_fit[[1]]$opt$aicc,ar_fit[[1]]$opt$aicc,asym_fit[[1]]$opt$aicc,curv_fit[[1]]$opt$aicc,X1_fit[[1]]$opt$aicc,X2_fit[[1]]$opt$aicc,body_fit[[1]]$opt$aicc,vol_fam_fit[[1]]$opt$aicc),
	OU = c(vol_fit[[2]]$opt$aicc,ar_fit[[2]]$opt$aicc,asym_fit[[2]]$opt$aicc,curv_fit[[2]]$opt$aicc,X1_fit[[2]]$opt$aicc,X2_fit[[2]]$opt$aicc,body_fit[[2]]$opt$aicc,vol_fam_fit[[2]]$opt$aicc),
	EB = c(vol_fit[[3]]$opt$aicc,ar_fit[[3]]$opt$aicc,asym_fit[[3]]$opt$aicc,curv_fit[[3]]$opt$aicc,X1_fit[[3]]$opt$aicc,X2_fit[[3]]$opt$aicc,body_fit[[3]]$opt$aicc,vol_fam_fit[[3]]$opt$aicc),
	WN = c(vol_fit[[4]]$opt$aicc,ar_fit[[4]]$opt$aicc,asym_fit[[4]]$opt$aicc,curv_fit[[4]]$opt$aicc,X1_fit[[4]]$opt$aicc,X2_fit[[4]]$opt$aicc,body_fit[[4]]$opt$aicc,vol_fam_fit[[4]]$opt$aicc),
	lambda = c(vol_lambda$lambda,ar_lambda$lambda,asym_lambda$lambda,curv_lambda$lambda,X1_lambda$lambda,X2_lambda$lambda,body_lambda$lambda,vol_fam_lambda$lambda),
	lambda_p = c(vol_lambda$P,ar_lambda$P,asym_lambda$P,curv_lambda$P,X1_lambda$P,X2_lambda$P,body_lambda$P,vol_fam_lambda$P))

### Print the results
sink("model_fitting_table.txt")
print(fit_table)
sink()

### Run arbutus for best fitting models
arb_vol_EB <- arbutus(vol_fit[[3]])
pdf("arbutus_vol_EB.pdf")
print(plot(arb_vol_EB))
dev.off()

arb_ar_EB <- arbutus(ar_fit[[3]])
pdf("arbutus_ar_EB.pdf")
print(plot(arb_ar_EB))
dev.off()

save.image(file="egg_analysis_model_fitting.RData")