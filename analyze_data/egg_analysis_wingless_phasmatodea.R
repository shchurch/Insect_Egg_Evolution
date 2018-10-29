### Created by SHC April 2018

### These commands can be used to ancestrally reconstruct shifts in wingless phasmatodea

source("analyze_data/egg_analysis_build_dataframe.R")
source("analyze_data/egg_analysis_read_trees.R")
source("analyze_data/egg_analysis_eco_functions.R")

### Read in the wingless phasmatodea table
ecology_table <- read.delim("analyze_data/ecology_table_wingless_phasmatodea.csv",header=T,stringsAsFactors=F)
eco_regimes <- c("NA","winged","wingless","partial_wing")
eco_data <- data.frame(name = ecology_table$name, rank = ecology_table$rank)
eco_data$ecology <- ecology_table$wing
default_regime <- "NA"
egg_database <- egg_database %>% filter(order == "Phasmatodea")
analysis_name <- "ecology_asr_wingless_Phasmatodea"

eco_data[which(eco_data$ecology == "partial_wing"),]$ecology <- "wingless"

### Set up the egg + ecology dataframe
# Randomly order the rows
egg_eco_data <- egg_database[sample(nrow(egg_database)),]
egg_eco_data$eco_regime <- rep(default_regime,nrow(egg_eco_data))

### Select a tree
egg_eco_data$rank <- egg_eco_data$genus
tree <- genus_mcc_tree

### Update database by ecology table
### Proceed from largest taxonomic category to smallest
egg_eco_data <- update_regimes_multistate("order",egg_eco_data$order,eco_regimes)
egg_eco_data <- update_regimes_multistate("suborder",egg_eco_data$suborder,eco_regimes)
egg_eco_data <- update_regimes_multistate("superfamily",egg_eco_data$superfamily,eco_regimes)
egg_eco_data <- update_regimes_multistate("family",egg_eco_data$family,eco_regimes)
egg_eco_data <- update_regimes_multistate("subfamily",egg_eco_data$subfamily,eco_regimes)
egg_eco_data <- update_regimes_multistate("tribe",egg_eco_data$tribe,eco_regimes)
egg_eco_data <- update_regimes_multistate("genus",egg_eco_data$genus,eco_regimes)
egg_eco_data <- update_regimes_multistate("species",egg_eco_data$name,eco_regimes)

### Apply filter for NA (used with Phasmatodea wingless dataset)
egg_eco_data <- egg_eco_data %>% filter(!(eco_regime == "NA"))

### make regimes numeric
egg_eco_data$discrete_regime <- as.numeric(as.factor(egg_eco_data$eco_regime))

