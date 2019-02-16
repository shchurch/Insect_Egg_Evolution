### Created by SCH April 2018

### These functions are used to execute phylogenetic regressions on traits of interest

library(nlme)

### This function runs a pgls analysis on a set of traits with a tree
pgls.trees <- function(data,tree,form){
  # The tree is pruned to include tips in the dataset
  pruned <- drop.tip(tree,setdiff(tree$tip.label,rownames(data)))
  # The dataset is trimmed to include only those in the tree
  red_dat <- data[!rownames(data) %in% setdiff(rownames(data),pruned$tip.label),]
  if(exists("blom_flag")) {
    # build a corBlomberg (ACDC / EB model) correlation structure  
    lm_str <- corBlomberg(value=1.3,phy=pruned,fixed=T)
  } else {
    # build a Brownian Motion correlation struction
    lm_str <- corBrownian(phy=pruned)
  }
  # run regression with model
  outskis <- gls(as.formula(form),correlation = lm_str,data = red_dat,method="ML")
  #outskis <- lm(form,dat)
  return(summary(outskis))
}

### This function builds a dataframe for running a pgls analysis
### It expects a data frame with rank, trait1, and trait2
build_pgls <- function(data,tree) {
  # The tree is pruned to include tips in the dataset
  pruned <- drop.tip(tree,setdiff(tree$tip.label,data$rank))
  # The dataset is trimmed to include only those in the tree
  red_dat <- subset(data,!(rank %in% setdiff(data$rank,pruned$tip.label)))

  # The dataframe for pgls is built        
  pgls_dat <- data.frame(trait1 = red_dat$trait1, trait2 = red_dat$trait2 ,rank = red_dat$rank)

  # Rows with NA and infinite values are removed
  is.na(pgls_dat) <- sapply(pgls_dat, is.infinite)
  pgls_dat <- na.omit(pgls_dat)

  # The dataset is grouped according to the taxonomic rank and a random entry is chosen
  pgls_dat <- data.frame(pgls_dat[sample(nrow(pgls_dat)),] %>% group_by(rank) %>% slice(1L))

  row.names(pgls_dat) <- pgls_dat$rank

  return(pgls_dat)
}

### This function builds a dataframe for running a pgls analysis against residuals
### It expects a data frame with rank, trait1, trait2, and indep (independent variable)
build_resid_pgls <- function(data,tree) {
  # The dataset is trimmed to include only those in the tree 
  red_dat <- subset(data,!(rank %in% setdiff(data$rank,tree$tip.label)))

  # The dataframe for the pgls is built        
  pgls_dat <- data.frame(trait1 = red_dat$trait1, trait2 = red_dat$trait2, indep = red_dat$indep, rank=red_dat$rank)

  # Rows with NA and infinite values are removed
  is.na(pgls_dat) <- sapply(pgls_dat,is.infinite)
  pgls_dat <-na.omit(pgls_dat)

  # The dataset is grouped according to the taxonomic rank and a random entry is chosen
  pgls_dat <- data.frame(pgls_dat[sample(nrow(pgls_dat)),] %>% group_by(rank) %>% slice(1L))
  
  # The dataframe for the pgls against the independent variable is built        
  Y <- data.frame(trait1 = pgls_dat$trait1, trait2 = pgls_dat$trait2, row.names = pgls_dat$rank)
  X <- pgls_dat$indep
  names(X) <- pgls_dat$rank
  
  # The tree is pruned to include tips in the dataset
  pruned <- drop.tip(tree,tip=setdiff(tree$tip.label,pgls_dat$rank))

  # The residuals are calculated with Brownian Motion
  resid <- phyl.resid(pruned,X,Y,method="BM")
  
  # The dataframe for the pgls of the residuals is built
  resid_dat <- data.frame(trait1 = resid$resid[,1],trait2 = resid$resid[,2],rank = rownames(resid$resid))

  return(resid_dat)
  # This function returns a dataframe for examining residuals
  #resid_plus_trait <- left_join(resid_dat,pgls_dat,by="rank") %>% rename(resid_trait1 = trait1.x, resid_trait2 = trait2.x, trait1 = trait1.y, trait2 = trait2.y)
}

### This function sets up and runs a pgls on all taxa and across predetermined groups
### It expects a dataframe with rank, trait1, and trait2
run_all_taxa_and_by_group_pgls <- function(data,tree,name,group_list){
  # Build the pgls dataset across all taxa
  all_dat <- build_pgls(data,tree)
  # Run pgls across all taxa
  all_pgls_out <- pgls.trees(all_dat,tree,"trait2 ~ trait1")

  # Create group datasets
  by_group_dat <- lapply(group_list,function(x) { data %>% filter(group == x) })

  # Build the pgls dataset across groups
  by_group_dat <- lapply(by_group_dat,build_pgls,tree = tree)
  # Run pgls across groups
  by_group_pgls_out <- lapply(by_group_dat,pgls.trees,tree = tree,form = "trait2 ~ trait1")

  # Add group names to results
  all_dat$group <- plyr::mapvalues(all_dat$rank,data$rank,data$group,warn_missing=F)

  return(list(all_pgls_out,by_group_pgls_out))
}

### This function sets up and runs a residual pgls on all taxa and across predetermined groups
### It expects a dataframe with rank, trait1, and trait2
run_all_taxa_and_by_group_resid_pgls <- function(data,tree,name,group_list){
  # Build the pgls residual dataset across all taxa
  all_dat <- build_resid_pgls(data,tree)
  # Runt the pgls residual across all taxa
  all_pgls_out <- pgls.trees(all_dat,tree,"trait2 ~ trait1")

  # Create group datasets
  by_group_dat <- lapply(group_list,function(x) { data %>% filter(group == x) })

  # Build the pgls dataset across groups
  by_group_dat <- lapply(by_group_dat,build_resid_pgls,tree = tree)
  # Run pgls across groups
  by_group_pgls_out <- lapply(by_group_dat,pgls.trees,tree = tree,form = "trait2 ~ trait1")

  # Add group names to results
  all_dat$group <- plyr::mapvalues(all_dat$rank,data$rank,data$group,warn_missing=F)

  return(list(all_pgls_out,by_group_pgls_out))
}

### This function sets up and runs a pgls on all taxa and across a list of regimes determined by group_list
### It expects a dataframe with rank, trait1, and trait2
run_eco_regime_pgls <- function(data,tree,name,group_list){
  # Build the pgls dataset across all taxa
  all_dat <- build_pgls(data,tree)
  # Run pgls across all taxa
  all_pgls_out <- pgls.trees(all_dat,tree,"trait2 ~ trait1")


  # Create group datasets
  by_group_dat <- lapply(group_list,function(x) { data %>% filter(group == x) })
  # Build the pgls dataset across groups
  by_group_dat <- lapply(by_group_dat,build_pgls,tree = tree)
  # Filter by groups with n > threshold
  threshold <- 20
  num_group <- data.frame(group = group_list, n = lapply(by_group_dat,nrow) %>% unlist())
  by_group_dat <- by_group_dat[ - which(num_group$n < threshold)]
  # Run pgls across groups
  by_group_pgls_out <- lapply(by_group_dat,pgls.trees,tree = tree,form = "trait2 ~ trait1")

  return(list(all_pgls_out,by_group_pgls_out))
}

### This function sets up and runs a pgls to plot regression lines
### It expects a dataframe with rank, trait1, and trait2
run_plotting_pgls <- function(data,tree,name){
  # Build the pgls dataset across all taxa
  all_dat <- build_pgls(data,tree)
  # Run pgls across all taxa
  all_pgls_out <- pgls.trees(all_dat,tree,"trait2 ~ trait1")

  # Create group datasets
  by_group_dat <- lapply(group_list,function(x) { data %>% filter(group == x) })

  # Build the pgls dataset across groups
  by_group_dat <- lapply(by_group_dat,build_pgls,tree = tree)
  # Run pgls across groups
  by_group_pgls_out <- lapply(by_group_dat,pgls.trees,tree = tree,form = "trait2 ~ trait1")

  # Add group names to results
  all_dat$group <- plyr::mapvalues(all_dat$rank,data$rank,data$group,warn_missing=F)

  # Plot results
  pdf(paste(name,".pdf",sep=""),width=6,height=6,useDingbats=F)
    g1 <- ggplot(all_dat,aes(x = trait1, y = trait2, color = group)) + geom_point() + scale_color_manual(values = mrk) + theme(legend.position = "none") #+ coord_fixed()
    g1 <- g1 + geom_abline(slope = 1,intercept = all_pgls_out$coefficients[1],size=1, linetype=2)
    #g1 <- g1 + geom_abline(slope = all_pgls_out$coefficients[2],intercept = all_pgls_out$coefficients[1],size=2)
    g1 <- g1 + geom_abline(slope = by_group_pgls_out[[1]]$coefficients[2],intercept = by_group_pgls_out[[1]]$coefficients[1],color = mrk[group_list[1]])
    g1 <- g1 + geom_abline(slope = by_group_pgls_out[[2]]$coefficients[2],intercept = by_group_pgls_out[[2]]$coefficients[1],color = mrk[group_list[2]])
    g1 <- g1 + geom_abline(slope = by_group_pgls_out[[3]]$coefficients[2],intercept = by_group_pgls_out[[3]]$coefficients[1],color = mrk[group_list[3]])
    g1 <- g1 + geom_abline(slope = by_group_pgls_out[[4]]$coefficients[2],intercept = by_group_pgls_out[[4]]$coefficients[1],color = mrk[group_list[4]])
    g1 <- g1 + geom_abline(slope = by_group_pgls_out[[5]]$coefficients[2],intercept = by_group_pgls_out[[5]]$coefficients[1],color = mrk[group_list[5]])
    g1 <- g1 + geom_abline(slope = by_group_pgls_out[[6]]$coefficients[2],intercept = by_group_pgls_out[[6]]$coefficients[1],color = mrk[group_list[6]])
    g1 <- g1 + geom_abline(slope = by_group_pgls_out[[7]]$coefficients[2],intercept = by_group_pgls_out[[7]]$coefficients[1],color = mrk[group_list[7]])
    print(g1)
  dev.off()
}

### This function sets up and runs a residual pgls to plot regression lines
### It expects a dataframe with rank, trait1, and trait2
run_plotting_resid_pgls <- function(data,tree,name){
  # Build the pgls residual dataset across all taxa
  all_dat <- build_resid_pgls(data,tree)
  # Runt the pgls residual across all taxa
  all_pgls_out <- pgls.trees(all_dat,tree,"trait2 ~ trait1")

  # Create group datasets
  by_group_dat <- lapply(group_list,function(x) { data %>% filter(group == x) })

  # Build the pgls dataset across groups
  by_group_dat <- lapply(by_group_dat,build_resid_pgls,tree = tree)
  # Run pgls across groups
  by_group_pgls_out <- lapply(by_group_dat,pgls.trees,tree = tree,form = "trait2 ~ trait1")

  # Add group names to results
  all_dat$group <- plyr::mapvalues(all_dat$rank,data$rank,data$group,warn_missing=F)

  # Plot the results
  pdf(paste(name,".pdf",sep=""),width=6,height=6,useDingbats=F)
    g1 <- ggplot(all_dat,aes(x = trait1, y = trait2, color = group)) + geom_point() + scale_color_manual(values = mrk) + theme(legend.position = "none")
    g1 <- g1 + geom_abline(slope = 1,intercept = all_pgls_out$coefficients[1],size=1, linetype=2)
    #g1 <- g1 + geom_abline(slope = all_pgls_out$coefficients[2],intercept = all_pgls_out$coefficients[1],size=2)
    g1 <- g1 + geom_abline(slope = by_group_pgls_out[[1]]$coefficients[2],intercept = by_group_pgls_out[[1]]$coefficients[1],color = mrk[group_list[1]])
    g1 <- g1 + geom_abline(slope = by_group_pgls_out[[2]]$coefficients[2],intercept = by_group_pgls_out[[2]]$coefficients[1],color = mrk[group_list[2]])
    g1 <- g1 + geom_abline(slope = by_group_pgls_out[[3]]$coefficients[2],intercept = by_group_pgls_out[[3]]$coefficients[1],color = mrk[group_list[3]])
    g1 <- g1 + geom_abline(slope = by_group_pgls_out[[4]]$coefficients[2],intercept = by_group_pgls_out[[4]]$coefficients[1],color = mrk[group_list[4]])
    g1 <- g1 + geom_abline(slope = by_group_pgls_out[[5]]$coefficients[2],intercept = by_group_pgls_out[[5]]$coefficients[1],color = mrk[group_list[5]])
    g1 <- g1 + geom_abline(slope = by_group_pgls_out[[6]]$coefficients[2],intercept = by_group_pgls_out[[6]]$coefficients[1],color = mrk[group_list[6]])
    g1 <- g1 + geom_abline(slope = by_group_pgls_out[[7]]$coefficients[2],intercept = by_group_pgls_out[[7]]$coefficients[1],color = mrk[group_list[7]])
    print(g1)
  dev.off()
}
