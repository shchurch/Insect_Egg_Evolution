### Created by SHC April 2018

### This code executes the allometry comparisons on the insect egg data
### This code generates the summary plots only
source("analyze_data/egg_analysis_build_dataframe.R")
source("analyze_data/egg_analysis_read_trees.R")
source("analyze_data/egg_analysis_pgls_functions.R")
library(phytools)

tree <- genus_mcc_tree

### set groups for subgroup analyses
group_list <- c("Hymenoptera","Condylognatha","Antliophora","Neuropteroidea","Amphiesmenoptera","Polyneoptera","Palaeoptera") 

### PGLS length vs width
length_width <- egg_database %>% select(genus,logX1,logX2,group) %>% rename(rank = genus, trait1 = logX1, trait2= logX2)
run_plotting_pgls(length_width,tree,"length_width")

### PGLS length vs asymmetry, width = independent
length_asym_width <- egg_database %>% select(genus,logX1,sqasym,logX2,group) %>% rename(rank = genus, trait1 = logX1, trait2= sqasym, indep = logX2)
run_plotting_resid_pgls(length_asym_width,tree,"length_asym_width")

### PGLS length vs curvature, width = independent
length_curv_width <- egg_database %>% select(genus,logX1,sqcurv,logX2,group) %>% rename(rank = genus, trait1 = logX1, trait2= sqcurv, indep = logX2)
run_plotting_resid_pgls(length_curv_width,tree,"length_curv_width")

### build body size dataset
# Set factor for resampling the body size dataframe, here no downsampling used
downsample_factor <- 1.0
fam_count_threshold <- 0
source("analyze_data/egg_analysis_body_size.R")

### PGLS egg size vs body size
vol_body <- egg_database_family_body %>% select(family,logvol,logbodyvol,group) %>% rename(rank = family, trait1 = logbodyvol, trait2= logvol)
run_plotting_pgls(vol_body,fam_tree,"vol_body")

### PGLS length vs width, body size = independent
length_width_body <- egg_database_family_body %>% select(family,logX1,logX2,logbody,group) %>% rename(rank = family, trait1 = logX1, trait2= logX2, indep = logbody)
run_plotting_resid_pgls(length_width_body,fam_tree,"length_width_body")

