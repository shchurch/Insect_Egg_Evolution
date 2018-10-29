### Created by SHC April 2018

### This code reads in the avian egg dataset from Stoddard et al 2017
### It builds comparable parameters between the insect and avian datasets

source("analyze_data/egg_analysis_build_dataframe.R")

### read in Stoddard et al Supplemental data
stoddard_supplement <- read.delim("analyze_data/stoddard.tsv",stringsAsFactors = F)

# this function calculates asymmetry as a ratio of quartile widths
# from a fitted lambda value
calc_q <- function(x) {
	val1 = (1 - x)
	val2 = (x - 1)
	val3 = (1 + x)
	q1 = val1 / val3
	q2 = val2 / val3
	q3 = (0.5)^q1
	q4 = (1.5)^q2
	q = q3 * q4
	return(q)
}

# convert bird parameters to bug parameters
bird_egg_data <- stoddard_supplement %>% 
		# rename variables according to egg database
		rename(X1 = AvgLength..cm.,
			el = Ellipticity,
			family = Family,
			order = Order,
			name = Species) %>% 
		# convert cm to mm for length
		mutate(X1 = X1 * 10,
			# calculate aspect ratio from ellipticity
			ar = el + 1,
			# calculate width
			X2 = X1 / ar,
			# calculate asymmetry as a ratio from lambda
			asym = calc_q(Asymmetry + 1) - 1,
			# set curvature to 0
			curv = 0,
			# calculate volume as an ellipsoid
			vol = ((4/3) * pi * (X1/2) * (X2/2)^2)) %>%
		# transform the data
		mutate(logar = log10(ar),
			logX1 = log10(X1),
			logX2 = log10(X2),
			logvol = log10(vol),
			sqcurv = 0,
			sqasym = sqrt(asym))

# Add a group label
bird_egg_data$group <- "Aves"

# Combine the avian and insect databases
insect_bird_egg_data <- rbind.fill(bird_egg_data, egg_database)
