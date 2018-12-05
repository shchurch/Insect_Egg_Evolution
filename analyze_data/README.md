# Evolutionary analysis in R

This file describes how to use the following R scripts to reproduce the figures in the Insect Egg Evolution manuscripts
All files in this directory were written by shchurch between 2017 and 2018

Each file should be executed in R. Package versions used in genearting these figures are listed at the end of this file

Some analyses rely on ecological and developmental data, taken from the primary literature. These data are provided as csv files, with sources listed as bibtex IDs, and full references provided in the file `ecology_development_size_references.bib`

## Create the cleaned egg database 

`egg_analysis_build_dataframe.R`

## Egg morphospace and model fitting 

```
egg_analysis_plot_morphospaces.R
egg_analysis_model_fitting.R
egg_analysis_rates_ancestral_states.R
```

### dependent files:

```
egg_analysis_read_trees.R
egg_analysis_body_size.R
egg_analysis_geiger_functions.R
egg_analysis_convert_stoddard.R
```

`groups_orders.txt` _(list of orders in each of 7 insect clades, e.g. Antliophora)_

`rainford_family_key.txt` _(list of taxonomic names in the Rainford et al. 2014 analysis)_

`family_body_sizes.csv` _(list of max and min body size for insect families, from Rainford et al. 2014)_

`stoddard.csv` _(avian egg size and shape data downloaded from Stoddard et al. 2014)_

`polyembryonic.csv` _(list of insect taxa described as polyembryonic with references)_

`bamm.config` _(configuration file to be updated with fitted parameters for bamm analyses)_

## Test allometric relationships 


`egg_analysis_allometry.R`
* rerun after recoding `corBlomberg` in `egg_analysis_pgls_functions.R`
* rerun after recoding Rainford trees in `egg_analysis_read_trees.R`

```
egg_analysis_plot_allometry.R
egg_analysis_simulate_allometry_data.R
egg_analysis_simulate_body.R
egg_analysis_clades_by_tips.R
```

### dependent files:

`egg_analysis_pgls_functions.R`


## Test relationship between development and eggs 

`egg_analysis_developement.R`
* rerun after recoding corBlomberg in `egg_analysis_pgls_functions.R`
* rerun after recoding Rainford trees in `egg_analysis_read_trees.R`

### dependent files:

`development.csv`


## Test relationship between ecology and eggs 

`egg_analysis_ouwie_ecology.R`
* run 4x for aquatic, parasitoid, wingless, and migratory
* rerun 2x parasitoid and aquatic after changing relaxed to strict  (`some` --> `ancestral`) in `egg_analysis_eco_functions.R`
* rerun aquatic version after setting `aquatic` to `combined aquatic + riparian` in `egg_analysis_aquatic.R`
* rerun parasitoid version after setting `parasitoid` to `combined endo + ecto`  in `egg_analysis_parasitoid.R`

`egg_analysis_plot_ecology.R`

`egg_analysis_separate_ecoallometry.R`
* run 2x for aquatic and parasitoid

### dependent files:

```
egg_analysis_eco_functions.R
egg_analysis_aquatic.R
egg_analysis_parasitoid.R
egg_analysis_wingless_phasmatodea.R
```

`egg_analysis_aquatic.csv`  _(list of insect taxa described as aquatic, riparian, semiaquatic with references)_

`egg_analysis_parasitoid.csv` _(list of insect taxa described as parasitoid with references)_

`egg_analysis_migratory_lepidoptera.csv` _(list of Lepidoptera taxa described as migratory with references)_

`egg_analysis_wingless_phasmatodea.csv` _(list of Phasmatodea taxa described as wingless / flightless with references)_



## Print tables and database statistics 

`egg_analysis_latex_tables.R` _(Directory containing results from other analysis is hardcoded)_

```
egg_analysis_art_stats.R
egg_analysis_database_statistics.R
```

### dependent files

`taxon_count.csv` _(estimated number of taxa in insect clades, from OTT)_


# General session info 

```
R version 3.4.2 (2017-09-28)
Platform: x86_64-apple-darwin15.6.0 (64-bit)
Running under: macOS Sierra 10.12.6

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base 
```

### sessionInfo():Egg morphospace and model fitting 

```
other attached packages:
[1] BAMMtools_2.1.6    arbutus_0.1        geiger_2.0.6       phytools_0.6-60   
[5] maps_3.3.0         phangorn_2.4.0     ape_5.2            bindrcpp_0.2.2    
[9] ggplot2_3.0.0      dplyr_0.7.6        plyr_1.8.4         RColorBrewer_1.1-2

loaded via a namespace (and not attached):
[1] deSolve_1.21            subplex_1.5-4           gtools_3.8.1           
[4] tidyselect_0.2.4        purrr_0.2.5             lattice_0.20-35        
[7] colorspace_1.3-2        expm_0.999-3            rlang_0.2.2            
[10] pillar_1.3.0            glue_1.3.0              withr_2.1.2            
[13] bindr_0.1.1             munsell_0.5.0           combinat_0.0-8         
[16] gtable_0.2.0            caTools_1.17.1.1        mvtnorm_1.0-8          
[19] coda_0.19-1             parallel_3.4.2          Rcpp_0.12.18           
[22] KernSmooth_2.23-15      scales_1.0.0            gdata_2.18.0           
[25] plotrix_3.7-3           clusterGeneration_1.3.4 scatterplot3d_0.3-41   
[28] gplots_3.0.1            fastmatch_1.1-0         mnormt_1.5-5           
[31] digest_0.6.17           animation_2.5           numDeriv_2016.8-1      
[34] grid_3.4.2              quadprog_1.5-5          bitops_1.0-6           
[37] magrittr_1.5            lazyeval_0.2.1          tibble_1.4.2           
[40] crayon_1.3.4            pkgconfig_2.0.2         MASS_7.3-50            
[43] Matrix_1.2-14           assertthat_0.2.0        reshape_0.8.7          
[46] R6_2.2.2                igraph_1.2.2            nlme_3.1-137           
[49] compiler_3.4.2   
```

### sessionInfo(): Egg allometry and development 

```
other attached packages:
[1] phylolm_2.6        geiger_2.0.6       nlme_3.1-137       phytools_0.6-60   
[5] maps_3.3.0         phangorn_2.4.0     ape_5.2            bindrcpp_0.2.2    
[9] ggplot2_3.0.0      dplyr_0.7.6        plyr_1.8.4         RColorBrewer_1.1-2

loaded via a namespace (and not attached):
[1] Rcpp_0.12.18            pillar_1.3.0            compiler_3.4.2         
[4] bindr_0.1.1             digest_0.6.17           tibble_1.4.2           
[7] gtable_0.2.0            lattice_0.20-35         pkgconfig_2.0.2        
[10] rlang_0.2.2             Matrix_1.2-14           fastmatch_1.1-0        
[13] igraph_1.2.2            parallel_3.4.2          mvtnorm_1.0-8          
[16] expm_0.999-3            coda_0.19-1             withr_2.1.2            
[19] globals_0.12.3          combinat_0.0-8          scatterplot3d_0.3-41   
[22] grid_3.4.2              tidyselect_0.2.4        deSolve_1.21           
[25] reshape_0.8.7           glue_1.3.0              listenv_0.7.0          
[28] R6_2.2.2                plotrix_3.7-3           future.apply_1.0.1     
[31] animation_2.5           purrr_0.2.5             magrittr_1.5           
[34] codetools_0.2-15        scales_1.0.0            MASS_7.3-50            
[37] assertthat_0.2.0        mnormt_1.5-5            future_1.9.0           
[40] colorspace_1.3-2        numDeriv_2016.8-1       quadprog_1.5-5         
[43] subplex_1.5-4           lazyeval_0.2.1          munsell_0.5.0          
[46] crayon_1.3.4            clusterGeneration_1.3.4
```

### sesssionInfo(): Allometry 

```
other attached packages:
[1] OUwie_1.50         lattice_0.20-35    corHMM_1.22        GenSA_1.1.7       
[5] nloptr_1.0.4       nlme_3.1-137       phytools_0.6-60    maps_3.3.0        
[9] phangorn_2.4.0     ape_5.2            bindrcpp_0.2.2     ggplot2_3.0.0     
[13] dplyr_0.7.6        plyr_1.8.4         RColorBrewer_1.1-2

loaded via a namespace (and not attached):
[1] gmp_0.5-13.2            Rcpp_0.12.18            pillar_1.3.0           
[4] compiler_3.4.2          bindr_0.1.1             tibble_1.4.2           
[7] gtable_0.2.0            pkgconfig_2.0.2         rlang_0.2.2            
[10] Matrix_1.2-14           fastmatch_1.1-0         igraph_1.2.2           
[13] parallel_3.4.2          expm_0.999-3            coda_0.19-1            
[16] Rmpfr_0.7-1             withr_2.1.2             nnet_7.3-12            
[19] combinat_0.0-8          scatterplot3d_0.3-41    grid_3.4.2             
[22] tidyselect_0.2.4        paleotree_3.1.0         reshape_0.8.7          
[25] glue_1.3.0              R6_2.2.2                plotrix_3.7-3          
[28] animation_2.5           corpcor_1.6.9           purrr_0.2.5            
[31] magrittr_1.5            scales_1.0.0            MASS_7.3-50            
[34] assertthat_0.2.0        mnormt_1.5-5            colorspace_1.3-2       
[37] numDeriv_2016.8-1       quadprog_1.5-5          lazyeval_0.2.1         
[40] munsell_0.5.0           crayon_1.3.4            clusterGeneration_1.3.4


sessionInfo(): latex tables and database stats 

other attached packages:
[1] ggplot2_3.0.0      dplyr_0.7.6        plyr_1.8.4         RColorBrewer_1.1-2
[5] gridExtra_2.3      xtable_1.8-3      

loaded via a namespace (and not attached):
[1] Rcpp_0.12.18     withr_2.1.2      crayon_1.3.4     assertthat_0.2.0
[5] grid_3.4.2       R6_2.2.2         gtable_0.2.0     magrittr_1.5    
[9] scales_1.0.0     pillar_1.3.0     rlang_0.2.2      lazyeval_0.2.1  
[13] bindrcpp_0.2.2   glue_1.3.0       munsell_0.5.0    purrr_0.2.5     
[17] compiler_3.4.2   colorspace_1.3-2 pkgconfig_2.0.2  tidyselect_0.2.4
[21] bindr_0.1.1      tibble_1.4.2  
```