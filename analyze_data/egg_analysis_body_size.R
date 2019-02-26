### Created by SHC April 2018

### This code builds the egg size + body size dataset
### Body size data for insect families is taken from the Rainford publication

### get family data from Rainford paper
family_body_size <- read.delim("analyze_data/family_body_sizes.tsv",stringsAsFactors = F)

# randomize the order of the egg database
egg_database_body_size <- egg_database[sample(nrow(egg_database)),]

# downsample the database to model the effects of the random averaging
### downsample_factor <- 0.50
### The downsample factor is going to be set by the parent R command, when this file is run from source()

# filter the data to entries with egg volume, group by family, and select a random sample
egg_database_body_size <- egg_database_body_size %>% filter(!(is.na(vol))) %>% group_by(family) %>% sample_frac(downsample_factor)

# reduce to only families with more than 1 egg observation
### fam_count_threshold <- 1
### The fam_count_threshold is going to be set by the parent R command, when this file is run from source()

fam_count <- egg_database_body_size %>% group_by(family) %>% summarize(n = n())
fam_count_n <- fam_count %>% filter(n > fam_count_threshold)
egg_database_body_size <- egg_database_body_size %>% filter(family %in% fam_count_n$family)

### log transform the dataframe
family_body_size <- family_body_size %>% mutate(
			logxt = log10(xt),
			logmt = log10(mt),
			body = (xt + mt) / 2,
			logbody = (logxt + logmt) / 2,
			bodyvol = (10^(logbody))^3,
			logbodyvol = log10(bodyvol))

### filter out the order level rows
fam <- family_body_size %>% filter(!(family == ""))

### summarize egg data by family
all_fam <- egg_database_body_size %>% filter(family %in% fam$family) %>%
		group_by(family) %>% 
		summarize_at(vars(vol,ar,logX1,logX2,logar,logvol,X1,X2,asym,curv,sqcurv,sqasym),funs(mean(.,na.rm=T)))

### create full dataset of egg and body size by family
egg_database_family_body <- merge(fam,all_fam,by="family")
egg_database_family_body <- egg_database_family_body %>% mutate(ind = vol / bodyvol, logind = log10(ind), logcub = log10((vol)^(1/3)))
egg_database_family_body$rank <- egg_database_family_body$family

### get order level rows
ordr <- family_body_size %>% filter(family == "") 

### summarize egg data by order for selected orders
all_ord <- egg_database_body_size %>% filter(order %in% ordr$order) %>% 
		group_by(order) %>% 
		summarize_at(vars(vol,ar,logX1,logX2,logar,logvol,X1,X2,asym,curv,sqcurv,sqasym),funs(mean(.,na.rm=T)))

### combine order level and family level datasets
full_ord <- merge(ordr,all_ord,by="order")
full_ord <- full_ord %>% mutate(ind = vol / bodyvol, logind = log10(ind), logcub = log10((vol)^(1/3)))
full_ord$rank <- full_ord$order

egg_database_family_body <- rbind(egg_database_family_body,full_ord)

### add in groups
egg_database_family_body$group <- plyr::mapvalues(egg_database_family_body$order,groups$order,groups$group,warn_missing=F)


