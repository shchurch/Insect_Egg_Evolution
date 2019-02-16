### Created by SHC, April 2018

### This code prints out the statistics on the insect egg database

source("analyze_data/egg_analysis_build_dataframe.R")
source("analyze_data/egg_analysis_read_trees.R")
library(gridExtra)
library(xtable)
options(scipen=1000000)

total_entries <- nrow(egg_database_raw)

#############################
### MAKE BIB FILE

unique_bibs <- egg_database_raw %>% filter(!(duplicated(b)))
all_bibs <- data.frame("bib_ID" = unique_bibs$b,"journal" = unique_bibs$bib_journal,"author" = unique_bibs$bib_author,"year" = unique_bibs$bib_year,"title" = unique_bibs$bib_title)
#write.table(all_bibs,file="all_references_in_database.csv",sep="\t",row.names=F)

#############################
### BIBLIOGRAPHIC STATISTICS
all_done_bib_ids <- read.delim("all_done_bib_ids.txt",header=F,stringsAsFactors=F)
num_refs <- length(unique(all_done_bib_ids$V1))
num_examined <- length(unique(egg_database_raw$b))
num_journals <- length(unique(egg_database_raw$bib_journal))
num_authors <- length(unique(egg_database_raw$bib_author))

bib_statistics <- data.frame("number" = c(num_refs,num_examined,num_authors,num_journals),row.names = c("references examined", "references with egg information", "unique authors", "unique journals / books"))

refs_by_year <- ggplot(egg_database %>% filter(!duplicated(bibtex)),aes(x=year,fill=group)) +  geom_histogram(binwidth=1) + scale_fill_manual(values=mrk) + scale_x_continuous(breaks = seq(from = 1600, to = 2017, by = 10)) + theme(legend.position = "none") + xlab("year of publication") + ylab("number of publications")
entries_by_year <- ggplot(egg_database,aes(x = year, fill = group)) +geom_histogram(binwidth=1) + scale_fill_manual(values=mrk) + scale_x_continuous(breaks = seq(from = 1600, to = 2017, by = 10)) + theme(legend.position = "none") + xlab("year of publication") + ylab("number of entries")

print(xtable(bib_statistics,digits=0),file="bibliographic_statistics_latex.txt")


#############################
### TAXONOMY STATISTICS

total_entries <- nrow(egg_database_raw)
entries_with_taxonomic_information <- nrow(egg_database)
problem_order <- nrow(egg_database_raw %>% filter(problem == "order"))
problem_noname <- nrow(egg_database_raw %>% filter(problem == "no_name"))
problem_notax <- nrow(egg_database_raw %>% filter(problem == "no_taxonomy"))
problem_nospecies <- nrow(egg_database_raw %>% filter(problem == "no_species"))
entries_with_genus_species <- nrow(egg_database %>% filter(!(species == "")))
entries_with_just_genus <- nrow(egg_database %>% filter(species == ""))

get_taxonomy_statistics <- data.frame("number" = c(total_entries,entries_with_taxonomic_information,problem_noname,problem_notax,problem_order,entries_with_just_genus,entries_with_genus_species,problem_nospecies),row.names=c("total entries","entries with taxonomic information","entries with unresolvable name","entries with no taxonomic information found","entries with incompatible order found","entries with only genus taxonomy","entries with genus and species taxonomy","entries where genus could be resolved, species could not"))

entries_with_genus_species <- egg_database %>% filter(!(species == ""))
taxonomic_statistics_by_order <- egg_database %>% group_by(order) %>% summarise(entries = length(name),unique_families = length(unique(family)),unique_genera = length(unique(genus)),unique_species = length(unique(name)))
taxonomic_statistics_by_group <- egg_database %>% group_by(group) %>% summarise(entries = length(name),unique_families = length(unique(family)),unique_genera = length(unique(genus)),unique_species = length(unique(name)))
taxonomic_statistics <- data.frame("number" = c(length(unique(entries_with_genus_species$name)),length(unique(egg_database$genus)),length(unique(egg_database$family)),length(unique(egg_database$order))),row.names=c("unique hexapod species","unique hexapod genera","unique hexapod families","unique hexapod orders"))

print(xtable(taxonomic_statistics,digits=0),file="taxonomic_statistics_latex.txt")

#############################
### NON ELLIPSOIDAL EGGS

# calculate the difference in volume as rotational ellipsoid or triaxial ellipsoid
# filter out those that already have volume calculated
egg_database_breadth <- egg_database %>% filter(vol != "") %>%
		# filter out those without breadth data
		filter(X3 != "") %>% mutate(
        # calculate width difference
        width_diff = X2 - X3,
        perc_width_diff = (width_diff / X2),
        # calculate volume 2 ways
        vol_width_breadth = vol,
        logvol_width_breadth = logvol,
        vol_width_width = ((4/3) * pi * (X1/2) * (X2/2)^2),
        logvol_width_width = log10(vol_width_width),
        # calculate volume difference
        vol_diff = abs(vol_width_breadth - vol_width_width),
        perc_vol_diff = (vol_diff / vol))

# build table
entries_with_breadth <- nrow(egg_database_breadth)
breadth_table <- data.frame("number" = c(total_entries,entries_with_breadth),row.names=c("total entries","entries with breadth recorded"))

# plot difference as point size in PCA plot
breadth_scatter <- ggplot(egg_database,aes(x = logar, y = logvol)) + geom_point(data = egg_database, color = "grey", size = 1) + geom_point(data = egg_database_breadth,aes(color = group,size = perc_width_diff)) + scale_color_manual(values = mrk) + scale_size_continuous(range = c(1,2.5)) + theme(legend.position = "None")

# plot histogram of percent difference in width
breadth_difference_statistics <- ggplot(egg_database_breadth %>% filter(perc_width_diff > 0.001),aes(perc_width_diff)) + geom_histogram(bins=30)  + scale_x_log10(limits=c(0.001,10),breaks=c(0,0.01,0.1,1.0,10),labels=scales::percent) + theme(legend.position = "None") + ylab("number of entries") + xlab("percent difference") + ggtitle("width vs breadth")

# plot histogram of percent difference in volume
breadth_vol_difference_statistics <- ggplot(egg_database_breadth %>% filter(perc_vol_diff > 0.00001),aes(perc_vol_diff)) + geom_histogram(bins=30)  + scale_x_log10(limits=c(0.001,10),breaks=c(0,0.01,0.1,1.0,10),labels=scales::percent) + theme(legend.position = "None") + ylab("number of entries") + xlab("percent difference in volume") + ggtitle("D. Triaxial  vs rotationally symmetric ellipsoid")

#############################
### VARIATION WITHIN TEXT ENTRIES

# calculate the difference from the biggest to smallest eggs
# percent difference = max - min / median 
#                    OR 2* deviation / average
vol_range <- egg_database_raw %>% mutate(
            bigX1 = ifelse(!is.na(al), 
                (al + dl), 
                xl),
            smallX1 = ifelse(!is.na(al), 
                (al - dl), 
                ml),
            bigX2 = ifelse(!is.na(aw), 
                (aw + dw), 
                xw),
            smallX2 = ifelse(!is.na(aw), 
                (aw - dw), 
                mw),
            vol1 = ((4/3) * pi * (bigX1/2) * (bigX2/2)^2),
            vol2 =  ((4/3) * pi * (smallX1/2) * (smallX2/2)^2),
            percent_diff = ((abs(bigX1 - smallX1) / ((bigX1 + smallX1)/2))))

# plot histogram of percent difference within entries
plot_reported_variation_in_length <- ggplot(vol_range,aes(x = percent_diff)) + geom_histogram(bins=30) + scale_x_log10(limits=c(0.001,10),breaks=c(0,0.01,0.1,1.0,10),labels=scales::percent) + xlab("percent difference length") + ggtitle("A. Reported intraspecific variation") + ylab("number of entries") + theme(legend.position = "none")

#############################
### VARIATION ACROSS INDEPENDENT OBSERVATIONS

#select entries which have text lengths
egg_X1 <- egg_database_w_duplicates %>% filter(!is.na(txtX1))
#get those names which have more than one entry
dup <- which(duplicated(egg_X1$name))
dup_names <- egg_X1[dup,]$name
#get those entries
duplicated_names <- egg_X1[which(egg_X1$name %in% dup_names),]
#filter the ones with genus + species names, and select the length
duplicated_species <- duplicated_names %>% filter(!(species == "")) %>% select(name,X1)

#calculate the difference in reported volumes ( max - min )
diffs <- duplicated_species %>% group_by(name) %>% summarise(diff = max(X1) - min(X1))
#find those which have no difference - likely repeated entries
repeats <- which(diffs$diff < (0.1 * 10^-5))
repeat_names <- diffs[repeats,]$name
#find those with reported difference - likely more sampling
multiple_samples <- which(diffs$diff > (0.1 * 10^-5))
multiple_names <- diffs[multiple_samples,]$name
#calcuate the percent difference in length for those which are multiply sampled
percent_diffs <- duplicated_species %>% filter(name %in% multiple_names) %>% 
                    group_by(name) %>%       
                    mutate(X11 = max(X1),X12 = min(X1)) %>%    
                    mutate(percent_diff = ((abs(X11 - X12) / ((X11 + X12)/2))))

dup_groups <- egg_database_w_duplicates %>% filter(name %in% percent_diffs$name) %>% select(name)
differences_in_reported_text_lengths <- percent_diffs %>% select(name,percent_diff) %>% left_join(dup_groups,by="name")

# plot the percent difference between observations
plot_diff_reported_txt_X1 <- ggplot(differences_in_reported_text_lengths,aes(percent_diff)) + geom_histogram(bins=30) + scale_x_log10(limits=c(0.001,10),breaks=c(0,0.01,0.1,1.0,10),labels=scales::percent) + xlab("percent difference in length") + ggtitle("B. Intraspecific variation across entries") + ylab("number of entries") + theme(legend.position = "none")

# write table
duplicated_statistics_table <- data.frame(entries = c(length(unique(entries_with_genus_species$name)),length(diffs$name),length(repeats),length(multiple_samples)),row.names = c("unique species in the database","species with multiple entries", "observations repeated across publications", "species with multiple independent observations"))

#############################
### VARIATION BETWEEN IMAGES AND TEXT

# calculate image volume and percent difference
ims_txt <- egg_database %>% 
                filter(!is.na(txtX1)) %>%
                filter(!is.na(imX1)) %>%
                select(name,txtX1,imX1,group,ID) %>% 
                mutate(X11 = txtX1, X12 = imX1) %>%
                mutate(percent_diff = ((abs(X11 - X12) / ((X11 + X12)/2))))

differences_in_lengths_text_images <- ims_txt %>% select(name,percent_diff,group)

# plot histogram of percent difference between text and images
plot_diff_txt_ims_X1 <- ggplot(differences_in_lengths_text_images,aes(percent_diff)) + geom_histogram(bins=30) + scale_x_log10(limits=c(0.001,10),breaks=c(0,0.01,0.1,1.0,10),labels=scales::percent) + xlab("percent difference in length") + ggtitle("C. Text description vs measured image") + ylab("number of entries") + theme(legend.position = "none")


#############################
### PRECISION IN TEXT ENTRIES

# flag those which are greater than 10 mm in length
bigger_than_10 <- which(egg_database_raw$pl > 10)
smaller_than_10 <- which(egg_database_raw$pl < 10)

# read in the data preserving trailing 0s
data_char <- read.csv("egg_database.csv",sep="\t",colClasses = "character")
# select just lengths
pl_char <- data_char$pl
# remove decimals
pl_char <- gsub("\\.","",pl_char)
# calculate number of significant figures (width of measurement in characters)
sigfigs <- nchar(pl_char)
# reduce those larger than 10 by 1
sigfigs[bigger_than_10] <- sigfigs[bigger_than_10] - 1

# translate the nubmer of significant figures into units
units <- plyr::mapvalues(sigfigs,c(1,2,3,4,5,6,7),c(1,0.1,0.01,0.001,0.0001,0.00001,0.000001))
# correct those larger than 1
units[bigger_than_10] <- units[bigger_than_10] * 0.1
# divide the observed length by the number of units used to measure
units_div <- units / egg_database_raw$pl
# calculate number of significant figures (width of measurement in characters)
ud <- nchar(units_div)

# change significant figures to decimal places
precision_absolute_units <- data.frame(sigfigs = sigfigs) %>% filter(sigfigs != 0) %>% mutate(decimals = sigfigs -1)
precision_relative_units <- data.frame(units_div = units_div) %>% na.omit()

plot_precision_absolute <- ggplot(precision_absolute_units,aes(x = decimals)) + 
    geom_histogram(binwidth=0.5) + 
    xlab("decimals") + 
    ggtitle("F. Precision of entries in milimeters") + 
    ylab("number of entries") + 
    theme(legend.position = "none") + 
    scale_x_continuous(breaks=c(0:10))

# set up new percent scale
NRpercent <- function(x) {
    paste0(sapply(x * 100, scales::comma), "%")
}


plot_precision_relative <- ggplot(precision_relative_units,aes(x = units_div)) + 
    geom_histogram(binwidth=0.5) + 
    xlab("percent of egg length") + 
    ggtitle("G. Relative precision of measurement units") + 
    ylab("number of entries") + 
    theme(legend.position = "none") + 
    scale_x_log10(breaks = c(0.00001,0.0001,0.001,0.01,0.1,1,10,100),label = NRpercent)

plch <- gsub("^0+","",pl_char)
benford <- as.integer(substring(plch,1,1))
donoughe <- as.integer(substring(pl_char,nchar(pl_char)))

benfords_law <- ggplot() + aes(benford)+ geom_histogram(binwidth=0.5) + ylab("frequency") + xlab("first digit of egg length record") + scale_x_continuous(breaks = c(1:9))


donoughes_law <- ggplot() + aes(donoughe)+ geom_histogram(binwidth=0.5) + ylab("frequency") + xlab("last digit of egg length record") + scale_x_continuous(breaks = c(1:9))



######################
##### IMAGE TEXT TYPE STATISTICS

image_status_table <- table(egg_database_raw$im_status)

image_type_table <- data.frame("images" = c(nrow(egg_database_raw %>% filter(im_q_image_type ==  "r")),nrow(egg_database_raw %>% filter(im_q_image_type ==  "t")),nrow(egg_database_raw %>% filter(im_q_image_type ==  "s")),nrow(egg_database_raw %>% filter(im_q_image_type ==  "d"))),row.names=c("micrograph with reflected light","micrograph with transmitted light","scanning electron micrograph","drawing"))

ims_only <- egg_database %>% filter(is.na(txtvol)) %>% select(imvol)
txt_only <- egg_database %>% filter(is.na(imvol)) %>% select(txtvol)
entries_with_just_image <- nrow(na.omit(ims_only))
entries_with_just_text <- nrow(na.omit(txt_only))
entries_with_both <- nrow(egg_database %>% filter(!is.na(txtvol)) %>% filter(!is.na(imvol)))


im_type_year <- data.frame(type = egg_database_raw$im_q_image_type, year = egg_database_raw$bib_year)

text_entries_average <- nrow(egg_database_raw %>% filter(!is.na(al)))
text_entries_range <- nrow(egg_database_raw %>% filter(!is.na(ml)))
text_entries_point_estimate <- nrow(egg_database_raw %>% filter(!is.na(pl)))

vol1 <- nrow(egg_database_raw %>% filter(!is.na(av)))
vol2 <- nrow(egg_database_raw %>% filter(!is.na(mv)))
vol3 <- nrow(egg_database_raw %>% filter(!is.na(pv)))

text_entries_volume_only <- vol1 + vol2 + vol3

measurable_images <- nrow(egg_database %>% filter(!is.na(imvol)))

text_statistics_table <- data.frame(entries = c(text_entries_average, text_entries_range, text_entries_point_estimate,text_entries_volume_only),row.names = c("average and st.dev. of length", "range of length", "single point estimate of length", "no length - volume only"))

image_text_table <- data.frame("entries" = c(total_entries,
    sum(text_entries_average,text_entries_range,text_entries_point_estimate),
    text_entries_average,
    text_entries_range,
    text_entries_point_estimate,
    text_entries_volume_only,
    nrow(egg_database_raw %>% filter(!(im_status %in% c("","none","missing")))),
    measurable_images,entries_with_both),
    row.names=c("Total entries in egg database",
    "Entries with text description of length and width",
    "Length reported as average and deviation",
    "Length reported as range",
    "Single length value reported",
    "Only volume reported",
    "Entries with an image",
    "Images re-measured",
    "Entries with both text and image measurements"))

print(xtable(image_text_table,digits=0),file="image_text_latex.txt")


######################
##### SAMPLING PLOTS

#### NEED TO UPDATE THIS FILE
taxon_count <- read.delim("analyze_data/taxon_count.csv",sep=",",stringsAsFactors = F)
tax_orders <- taxon_count[which(taxon_count$OTL_rank == "order"),]

tax_orders <- data.frame("order" = tax_orders$taxon, "entries" = tax_orders$N_in_DB, "total_taxonomic_tips" = tax_orders$Ntips_tax) %>% 
        filter(order %in% groups$order)

taxord1 <- ggplot(tax_orders,aes(x=total_taxonomic_tips,y=entries)) + 
    geom_abline(linetype=2, size = 1.5, intercept=-2,slope=1,color="black") + 
    geom_point(size=2)+ 
    scale_x_log10(breaks=c(0,1,10,100,1000,10000,100000,1000000)) +
    scale_y_log10(breaks=c(0,1,10,100,1000,10000,100000,1000000)) +
    ggtitle("B") + 
    theme(legend.position="none",
    panel.grid.major = element_blank(),
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank()) + xlab("species per order") + ylab("entries per order")

tax_families <- taxon_count[which(taxon_count$taxon %in% egg_database_raw$tax_family),]
tax_families <- data.frame("family" = tax_families$taxon, "entries" = tax_families$N_in_DB, "total_taxonomic_tips" = tax_families$Ntips_tax) %>% 
        filter(family %in% egg_database$family)

taxfam1 <- ggplot(tax_families,aes(x=total_taxonomic_tips,y=entries))  + 
    geom_abline(linetype=2, size = 1.5, intercept=-2,slope=1,color="black") +
    geom_point(size=2,aes(text=family)) + 
    scale_x_log10(breaks=c(0,1,10,100,1000,10000,100000,1000000)) +
    scale_y_log10(breaks=c(0,1,10,100,1000,10000,100000,1000000)) +
    ggtitle("A") + 
    theme(legend.position="none",
    panel.grid.major = element_blank(),
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank()) + xlab("species per family") + ylab("entries per family")

pdf("Figure_sampling_statistics.pdf",height=5,width=15)
grid.arrange(taxfam1,taxord1,ncol=2)
dev.off()

#############################
### FIGURES FOR PUBLICATION

pdf("Figure_data_validation.pdf",height=10,width=15)
grid.arrange(plot_reported_variation_in_length,plot_diff_reported_txt_X1,plot_diff_txt_ims_X1,breadth_vol_difference_statistics,plot_precision_absolute,plot_precision_relative,ncol=2)
dev.off()