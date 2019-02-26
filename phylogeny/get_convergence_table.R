
library(tidyverse)
library(mcmcse)

get_ess = function(df){
  backbone = df$backbone
  order = df$order
  message(paste('Working on',backbone,order))
  parfiles = Sys.glob(str_c('mrbayes_smalltrees/treelinks_',backbone,'/',order,'*')) %>%
    Sys.readlink() %>%
    str_replace('.t$','.p') %>%
    str_replace('^','mrbayes_smalltrees/treelinks_misof/')
  names(parfiles) = 1:length(parfiles)
  
  plyr::ldply(parfiles, readr::read_tsv,skip=1,.id = 'run') %>%
    group_by(run) %>%
    mutate(index = 1:n()) %>%
    filter(index > quantile(index,0.1)) %>%
    ungroup %>%
    summarise(post_burnin_samples = n(),
              ESS_likelihood = ess(LnL),
              ESS_prior = ess(LnPr)) %>%
    mutate(backbone = backbone, order = order) %>%
    select(backbone, order, post_burnin_samples, ESS_likelihood, ESS_prior)
}

backbones = c('misof', 'rainford') 
orders = c('Archaeognatha',
          'Coleoptera',
          'Collembola',
          'Dermaptera',
          'Dictyoptera',
          'Diplura',
          'Diptera',
          'Ephemeroptera',
          'Grylloblattodea',
          'Hemiptera',
          'Hymenoptera',
          'Lepidoptera',
          'Mantodea',
          'Mantophasmatodea',
          'Mecoptera',
          'Megaloptera',
          'Neuroptera',
          'Odonata',
          'Orthoptera',
          'Phasmatodea',
          'Plecoptera',
          'Protura',
          'Psocodea',
          'Siphonaptera',
          'Strepsiptera',
          'Thysanoptera',
          'Trichoptera',
          'Zygentoma')

final_df = expand.grid(backbone = backbones, order = orders)

convergence_table = final_df %>%
  plyr::adply(1,get_ess)

readr::write_csv(convergence_table,'final_trees/convergence_table.csv')

