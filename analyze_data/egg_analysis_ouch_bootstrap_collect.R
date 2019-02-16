### This file was created by B. Medeiros starting Feb 2019

### We ran bootstraps for OUCH in parallel and will now collect results and put in a table
### The folder with bootstraps should be given as command-line arg
args = commandArgs(trailingOnly=TRUE)
infolder = file.path(args[1])

library(ouch)
library(tidyverse)
library(scales)
# we will not load plyr to avoid problems, but needs package plyr

### Function to read one bootstrap and return as data.frame
read_boot = function(inpath){
  load(inpath) #the object loaded is named bootrep
  return(bootrep) 
}

#now let's read all of them
files = list.files(infolder,full.names = T)
names(files) = 1:length(files)

boot_table = plyr::ldply(.data = files,
            .fun = read_boot,
            .id = 'replicate')

#and print confidence intervals nicely

#function to get confidence interval as text
get_ci = function(vec,p=0.05){
  vec %>%
    quantile(c(p/2,1-p/2)) %>%
    scales::number(accuracy = 0.0001) %>%
    str_c(collapse = ' -- ')
}

ci_table = boot_table %>%
  select(starts_with('optima'),
         starts_with('sigma'),
         starts_with('alpha')) %>%
  summarise_all(get_ci)

print('SIGMAS')
sigmas = ci_table %>%
  select(starts_with('sigma')) %>%
  unlist %>%
  matrix(nrow=4)
print(sigmas)

cat('\n\n\n')

print('ALPHAS')
alphas = ci_table %>%
  select(starts_with('alpha')) %>%
  unlist %>%
  matrix(nrow=4)
print(alphas)

cat('\n\n\n')

print('OPTIMA')
optima = ci_table %>%
  select(starts_with('optima')) %>%
  t %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  transmute(dependent_trait = rowname %>% 
              str_extract('trait[0-9]') %>%
              str_remove('[a-z]+'),
            independent_trait =  rowname %>% 
              str_extract('[0-9]$') %>%
              str_c('optimum_',.),
            value = V1) %>%
  spread(independent_trait,value)
print(optima)

outfile= file.path(infolder,'all_bootstraps.R')

save(file = outfile,boot_table,ci_table,sigmas,alphas,optima)

