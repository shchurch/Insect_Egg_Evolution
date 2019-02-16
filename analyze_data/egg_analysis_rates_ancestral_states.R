### Created by SHC May 2018

### This code reconstructs the ancestral state for egg traits
### It also sets up the parameters and plots for a BAMM analysis of evolutionary rates 

library(BAMMtools)
library(phytools)
source("analyze_data/egg_analysis_build_dataframe.R")
source("analyze_data/egg_analysis_read_trees.R")
source("analyze_data/egg_analysis_geiger_functions.R")

### select a tree
tree <- genus_mcc_tree

### Select a continuous trait to analyze
trait_of_interest <- "logvol"
#trait_of_interest <- "logar"
#trait_of_interest <- "sqasym"
#trait_of_interest <- "sqcurv"

### name the analysis
name <- paste("asr_rate",trait_of_interest,sep="_")

### Set up the data frame for analysis
dat <- egg_database[sample(nrow(egg_database)),] %>% 
	mutate(rank = genus,trait = eval(parse(text = trait_of_interest))) %>% 
	filter(!(is.na(trait))) %>% 
	group_by(rank) %>% slice(1L)

### BAMM

### Format the data frame for BAMM
bmm_red <- dat %>% mutate(taxa = rank) %>% 
		filter(rank %in% tree$tip.label) %>% 
		select(rank,trait)

### Use only taxa with phenotypic data
bmm <- drop.tip(tree,setdiff(tree$tip.label,bmm_red$rank))

### Write the files for BAMM
write.table(bmm_red,file="bamm_trait.txt",sep="\t",row.names=F,col.names=F,quote=F)
# Select appropriate priors (remember to change expected numbers of shifts to 10)
setBAMMpriors(phy = bmm, traits = "bamm_trait.txt")
write.tree(bmm,file="bamm_tree.tre")

#### OUTSIDE OF R, RUN BAMM
#### copy parameters over to bamm.config
#### set expected number of shifts to 10
#### bamm -c bamm.config

### Read in the BAMM results and plot the shifts on the tree
tree <- read.tree(file="bamm_tree.tre")
edata <- getEventData(tree, eventdata = paste("bamm_event_data.txt",sep=""), burnin=0.1, type="trait")
pdf(file=paste(name,"_ratesplot.pdf",sep=""),width=16,height=16)
	plot.bammdata(edata,lwd=2,legend=T,labels=T,cex=0.2,logcolor=T,breaksmethod="jenks",pal=c("#0D0887", "#CC4678", "#F0F921"))
dev.off()

### CONTINUOUS TRAIT ASR

### Format the data frame for phytools
ms <- build_geiger_dataset(dat = dat, tree = tree)
ms <- ms[order(match(names(ms),tree$tip.label))]

### Use only taxa with phenotypic data
ms_trim <- drop.tip(tree,setdiff(tree$tip.label,names(ms)))
ms_trim <- di2multi(ms_trim)

### Plot the reconstructed ancestral state
pdf(file=paste(name,".pdf",sep=""),height=16,width=16)
  # Perform the ancestral state reconstruction
  obj <- contMap(ms_trim,ms,plot=F)
  obj <- setMap(obj,colors=c("#440154","#21908C","#FDE725"))
  plot(obj,fsize=0.2,lwd=3,outline=F)
dev.off()

