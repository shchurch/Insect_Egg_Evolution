### Created by SHC April 2018

### These commands can be used to ancestrally reconstruct shifts in parasitism

source("analyze_data/egg_analysis_build_dataframe.R")
source("analyze_data/egg_analysis_read_trees.R")
source("analyze_data/egg_analysis_eco_functions.R")

### Read in the ecological table
ecology_table <- read.delim("analyze_data/ecology_table_aquatic.csv",header=T,stringsAsFactors=F)
eco_data <- data.frame(name = ecology_table$name, rank = ecology_table$rank)
default_regime <- "ancestral"

eco_data$ecology <- ecology_table$in_water
analysis_name <- "ecology_asr_in_water"

# Uncomment these to use aquatic + riparian, instead of aquatic oviposition
#eco_data$ecology <- ecology_table$combined
#analysis_name <- "ecology_asr_combined_aquatic"


### Set up the egg + ecology dataframe
# Randomly order the rows
egg_eco_data <- egg_database[sample(nrow(egg_database)),]
egg_eco_data$eco_regime <- rep(default_regime,nrow(egg_eco_data))

### Select a tree
egg_eco_data$rank <- egg_eco_data$genus
tree <- genus_mcc_tree

### Update database by ecology table
### Proceed from largest taxonomic category to smallest
egg_eco_data <- update_regimes_TF("order",egg_eco_data$order)
egg_eco_data <- update_regimes_TF("suborder",egg_eco_data$suborder)
egg_eco_data <- update_regimes_TF("superfamily",egg_eco_data$superfamily)
egg_eco_data <- update_regimes_TF("family",egg_eco_data$family)
egg_eco_data <- update_regimes_TF("subfamily",egg_eco_data$subfamily)
egg_eco_data <- update_regimes_TF("tribe",egg_eco_data$tribe)
egg_eco_data <- update_regimes_TF("genus",egg_eco_data$genus)
egg_eco_data <- update_regimes_TF("species",egg_eco_data$name)

### make regimes numeric
egg_eco_data$discrete_regime <- as.numeric(as.factor(egg_eco_data$eco_regime))
