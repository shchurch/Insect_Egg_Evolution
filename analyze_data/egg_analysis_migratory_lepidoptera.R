### Created by SHC April 2018

### These commands can be used to ancestrally reconstruct shifts in migratory lepidoptera

#source("analyze_data/egg_analysis_build_dataframe.R")
#source("analyze_data/egg_analysis_read_trees.R")
source("analyze_data/egg_analysis_eco_functions.R")

### Read in the wingless phasmatodea table
ecology_table <- read.delim("analyze_data/ecology_table_migratory.tsv",header=T,stringsAsFactors=F)
eco_regimes <- c("nonmigratory","migratory")
eco_data <- data.frame(name = ecology_table$name, rank = ecology_table$rank)
eco_data$ecology <- ecology_table$migratory
default_regime <- "ancestral"
analysis_name <- paste(analysis_name,"migratory_Lepidoptera",sep="_")
egg_database <- egg_database %>% filter(order == "Lepidoptera")

### Choose colors
eco_cols <- c("dark gray","#2056CE")

### Set up the egg + ecology dataframe
# Randomly order the rows
egg_eco_data <- egg_database[sample(nrow(egg_database)),]
egg_eco_data$eco_regime <- rep(default_regime,nrow(egg_eco_data))

### Select a taxonomic level
egg_eco_data$rank <- egg_eco_data$genus

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

