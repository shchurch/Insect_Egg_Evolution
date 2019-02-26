### Created by SHC  May 2018

### This script codes in species which are reported to lay eggs in ootheca
### It is based on egg_analysis_aquatic

source("analyze_data/egg_analysis_eco_functions.R")
### Read in the ecological table
ecology_table <- read.delim("analyze_data/ootheca.tsv",header=T,stringsAsFactors=F)
eco_data <- data.frame(name = ecology_table$name, rank = ecology_table$rank)
default_regime <- "ancestral"
eco_data$ecology <- ecology_table$ootheca

### There are no ambiguous taxa in the ootheca dataset
class_flag <- "relaxed"

egg_eco_data <- egg_database[sample(nrow(egg_database)),]
egg_eco_data$eco_regime <- rep(default_regime,nrow(egg_eco_data))
egg_eco_data$rank <- egg_eco_data$genus

egg_eco_data <- update_regimes_TF("order",egg_eco_data$order)
egg_eco_data <- update_regimes_TF("suborder",egg_eco_data$suborder)
egg_eco_data <- update_regimes_TF("superfamily",egg_eco_data$superfamily)
egg_eco_data <- update_regimes_TF("family",egg_eco_data$family)
egg_eco_data <- update_regimes_TF("subfamily",egg_eco_data$subfamily)
egg_eco_data <- update_regimes_TF("tribe",egg_eco_data$tribe)
egg_eco_data <- update_regimes_TF("genus",egg_eco_data$genus)
egg_eco_data <- update_regimes_TF("species",egg_eco_data$name)


