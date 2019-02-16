### This file was created by B. Medeiros starting Feb 2019

### This code reads an OUCH model and performs one replicate of parametric bootstrapping
### After running replicates in parallel, results will be read in another script
library(ouch)
library(stringr)

### This analysis will read SLURM array # and use it as random seed
cat(stringr::str_c('Random seed: ',Sys.getenv('SLURM_ARRAY_TASK_ID')))

### The R object with ouch model shoudl be given as command-line argument
### Let's load the R object
### We will keep in a separate environment to avoid name conflicts
args = commandArgs(trailingOnly=TRUE)
ouch_res = new.env()
load(args[1],envir = ouch_res)

### We will use input file name to name output files
prefix = basename(args[1])


### Now, let's do one bootstrap replicate
rseed = Sys.getenv('SLURM_ARRAY_TASK_ID')
bootrep = bootstrap(ouch_res$oum_fit,nboot = 1,seed = rseed)

### and save it
dir.create(path = file.path('ouch_boot', prefix),showWarnings = F,recursive = T)

outfile = file.path('ouch_boot',
                    prefix,
                    stringr::str_c('b_',Sys.getenv('SLURM_ARRAY_TASK_ID'),'.Rdata'))
save(bootrep,file = outfile)
