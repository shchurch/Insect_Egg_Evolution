source("analyze_data/egg_analysis_geiger_functions.R")

# read in arguments from the OUwie analysis wrapper
# 3 = testing 'parasitoid', 'internal', 'aquatic', 'combined_aquatic'

if(ecology_arg == "internal") {
	source("analyze_data/egg_analysis_parasitoid.R")

} else if(ecology_arg == "aquatic") {
	source("analyze_data/egg_analysis_aquatic.R")

} else if(ecology_arg == "parasitoid") {
	source("analyze_data/egg_analysis_parasitoid.R")

} else if(ecology_arg == "combined_aquatic") {
	source("analyze_data/egg_analysis_aquatic.R")

} else {
	warning("incorrect ecology argument")
	quit(status=1)
}

g <- build_geiger_dataset(egg_eco_data %>% mutate(rank = genus,trait = discrete_regime),tree)
reg_fit <- fit_discrete_model(g,tree,"ARD")
reg_model <- c(-reg_fit$opt$q12,reg_fit$opt$q12,reg_fit$opt$q21,-reg_fit$opt$q21)

simulated_discrete_trait <- rTraitDisc(phy=tree,model=matrix(reg_model,2),root.value="1",states=c("1","2"))

simulated_egg_eco_data <- data.frame(genus = names(simulated_discrete_trait),sim_discrete_regime = unname(simulated_discrete_trait),stringsAsFactors=F)

egg_eco_data <- inner_join(egg_eco_data,simulated_egg_eco_data,by="genus") %>% mutate(discrete_regime = sim_discrete_regime) %>% select(-eco_regime,-sim_discrete_regime)