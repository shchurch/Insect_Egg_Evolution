source("analyze_data/egg_analysis_build_dataframe.R")
source("analyze_data/egg_analysis_read_trees.R")

analysis_name <- "genome_size_egg_size"
args = commandArgs(trailingOnly=TRUE)
# This program takes two arguments
# 1 - which backbone tree to use, 'misof' [default] or 'rainford'
# 2 - which correlation matrix to use 'brownian' [default] or 'blomberg'
if(args[1] == "rainford") {
	genus_trees <- rainford_genus_trees
	genus_mcc_tree <- rainford_genus_mcc_tree
	analysis_name <- paste("rainford",analysis_name,sep="_")
}
if(args[2] == "blomberg") {
	blom_flag <- TRUE
	analysis_name <- paste("corBlomberg",analysis_name,sep="_")
}
source("analyze_data/egg_analysis_pgls_functions.R")

group_list <- c("Hymenoptera","Condylognatha","Antliophora","Neuropteroidea","Amphiesmenoptera","Polyneoptera","Palaeoptera")

genome_size <- read.delim("analyze_data/genome_size.tsv",header=T,stringsAsFactors=F)

genome_size$name <- paste(genome_size$genus,genome_size$species,sep="_")

run_genome_pgls <- function(tree) {
	genome <- genome_size %>% select(name,c_value) %>% mutate(logcval = log10(c_value))
	genome_vol <- left_join(egg_database,genome,by="name") %>% 
		filter(!(is.na(c_value))) %>% 
		rename(rank = genus, trait1 = logvol, trait2 = logcval)
	pgls_genome_vol <- run_all_taxa_and_by_group_pgls(genome_vol,tree,analysis_name,group_list)
	
	genome_genus <- genome_size[sample(nrow(genome_size)),] %>%
		select(genus,c_value) %>% 
		mutate(logcval = log10(c_value)) %>% 
		group_by(genus) %>% slice(1L)
	genome_vol_genus <- left_join(egg_database,genome_genus,by="genus") %>% 
		filter(!(is.na(c_value))) %>% 
		rename(rank = genus, trait1 = logvol, trait2 = logcval)
	pgls_genome_vol_genus <- run_all_taxa_and_by_group_pgls(genome_vol_genus,tree,paste(analysis_name,"genus",sep="_"),group_list)

	genome_vol_no_Polyneoptera <- left_join(egg_database %>% filter(group != "Polyneoptera"),genome,by="name") %>% 
		filter(!(is.na(c_value))) %>% 
		rename(rank = genus, trait1 = logvol, trait2 = logcval)
	pgls_genome_vol_no_Polyneoptera <- run_all_taxa_and_by_group_pgls(genome_vol_no_Polyneoptera,tree,paste(analysis_name,"without_Polyneoptera",sep="_"),
		c("Hymenoptera","Condylognatha","Antliophora","Neuropteroidea","Amphiesmenoptera","Palaeoptera"))

	return(list(pgls_genome_vol,pgls_genome_vol_genus,pgls_genome_vol_no_Polyneoptera))
}


genome_pgls_distribution_raw <- lapply(genus_trees,run_genome_pgls)

### Build a results table for all insects
get_genome_pgls_all_taxa_table <- function(pgls) {
	table <- round(data.frame(
	slope_min = min(sapply(genome_pgls_distribution_raw,function(x) {x[[pgls]][[1]]$coefficients[[2]]})),
	slope_max = max(sapply(genome_pgls_distribution_raw,function(x) {x[[pgls]][[1]]$coefficients[[2]]})),
	int_min = min(sapply(genome_pgls_distribution_raw,function(x) {x[[pgls]][[1]]$coefficients[[1]]})),
	int_max = max(sapply(genome_pgls_distribution_raw,function(x) {x[[pgls]][[1]]$coefficients[[1]]})),
	pval_min = min(sapply(genome_pgls_distribution_raw,function(x) {x[[pgls]][[1]]$tTable[[8]]})),
	pval_max = max(sapply(genome_pgls_distribution_raw,function(x) {x[[pgls]][[1]]$tTable[[8]]})),
	sample_size = min(sapply(genome_pgls_distribution_raw,function(x) {x[[pgls]][[1]]$dims$N}))),
	4)
	return(table)
}

genome_vol_table <- get_genome_pgls_all_taxa_table(1)
genome_vol_genus_table <- get_genome_pgls_all_taxa_table(2)
genome_vol_no_Polyneoptera_table <- get_genome_pgls_all_taxa_table(3)

### Build a results table for each group
get_group_genome_pgls_table <- function(pgls,value) {
	table <- round(data.frame(
	slope_min = min(sapply(genome_pgls_distribution_raw,function(x) {x[[pgls]][[2]][[value]]$coefficients[[2]]})),
	slope_max = max(sapply(genome_pgls_distribution_raw,function(x) {x[[pgls]][[2]][[value]]$coefficients[[2]]})),
	int_min = min(sapply(genome_pgls_distribution_raw,function(x) {x[[pgls]][[2]][[value]]$coefficients[[1]]})),
	int_max = max(sapply(genome_pgls_distribution_raw,function(x) {x[[pgls]][[2]][[value]]$coefficients[[1]]})),
	pval_min = min(sapply(genome_pgls_distribution_raw,function(x) {x[[pgls]][[2]][[value]]$tTable[[8]]})),
	pval_max = max(sapply(genome_pgls_distribution_raw,function(x) {x[[pgls]][[2]][[value]]$tTable[[8]]})),
	sample_size = min(sapply(genome_pgls_distribution_raw,function(x) {x[[pgls]][[2]][[value]]$dims$N}))),
	4)
	return(table)
}

genome_vol_group_table <- lapply(seq(length(group_list)),get_group_genome_pgls_table,pgls=1)
genome_vol_genus_group_table <- lapply(seq(length(group_list)),get_group_genome_pgls_table,pgls=2)

### Print the results
sink(paste(analysis_name,".txt",sep=""))
cat("Across all insects\n")
cat("\ngenome_size_egg_size\n")
print(genome_vol_table)
cat("\ngenome_size_egg_size_genus\n")
print(genome_vol_genus_table)
cat("\ngenome_size_egg_size_without_Polyneoptera\n")
print(genome_vol_no_Polyneoptera_table)
for(i in seq(1:length(group_list))) {
	cat("\n")
	print(group_list[i])
	cat("\ngenome_size_egg_size\n")
	print(genome_vol_group_table[[i]])
	cat("ngenome_size_egg_size_genus\n")
	print(genome_vol_genus_group_table[[i]])
}
sink()


tree <- genus_mcc_tree

genome <- genome_size %>% select(name,c_value) %>% mutate(logcval = log10(c_value))
genome_vol <- left_join(egg_database,genome,by="name") %>% 
	filter(!(is.na(c_value))) %>% 
	rename(rank = genus, trait1 = logvol, trait2 = logcval)
run_plotting_pgls(genome_vol,tree,paste(analysis_name,"plotting_pgls",sep="_"))

pdf(paste(analysis_name,".pdf",sep=""),width=6,height=6)
print(ggplot(genome_vol %>% select(rank,c_value,vol,group) %>% na.omit() %>% group_by(rank) %>% slice(1L),aes(x = vol,y = c_value,color = group)) + 
	geom_point() + 
	scale_x_log10() + 
	scale_y_log10() + 
	theme(legend.position = "none") + 
	scale_color_manual(values = mrk)) + 
	xlab("egg volume") + 
	ylab("genome C value")
dev.off()

genome_genus <- genome_size[sample(nrow(genome_size)),] %>%
		select(genus,c_value) %>% 
		mutate(logcval = log10(c_value)) %>% 
		group_by(genus) %>% slice(1L)
genome_vol_genus <- left_join(egg_database,genome_genus,by="genus") %>% 
		filter(!(is.na(c_value))) %>% 
		rename(rank = genus, trait1 = logvol, trait2 = logcval)
run_plotting_pgls(genome_vol_genus,tree,paste(analysis_name,"plotting_pgls",sep="_"))

pdf(paste(analysis_name,"_genus",".pdf",sep=""),width=6,height=6)
print(ggplot(genome_vol_genus %>% select(rank,c_value,vol,group) %>% na.omit() %>% group_by(rank) %>% slice(1L),aes(x = vol,y = c_value,color = group)) + 
	geom_point() + 
	scale_x_log10() + 
	scale_y_log10() + 
	theme(legend.position = "none") + 
	scale_color_manual(values = mrk)) +
	xlab("egg volume") + 
	ylab("genome C value")
dev.off()

pdf(paste(analysis_name,"_without_Polyneoptera",".pdf",sep=""),width=6,height=6)
print(ggplot(genome_vol %>% filter(group != "Polyneoptera") %>% select(rank,c_value,vol,group) %>% na.omit() %>% group_by(rank) %>% slice(1L),aes(x = vol,y = c_value,color = group)) + 
	geom_point() + 
	scale_x_log10() + 
	scale_y_log10() + 
	theme(legend.position = "none") + 
	scale_color_manual(values = mrk)) +
	xlab("egg volume") + 
	ylab("genome C value")
dev.off()

save.image(file=paste("egg_analysis",analysis_name,"workspace.RData",sep="_"))


