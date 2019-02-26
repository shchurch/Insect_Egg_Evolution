### Created by SHC, April 2018

### This code analyzes the results of the image parsing sensitivity test

library(RColorBrewer)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(xtable)
theme_set(theme_classic(base_size=18))

# Read in the results of the reclicking assay
art_egg_database_raw <- read.delim("analyze_data/art_eggs_database.tsv",stringsAsFactors=T)

# Filter out hte results in which aspect ratio was equal to or below 1 with curvature
# Axiomatically curvature cannot be deteceted in a sphere
# We additional define curvature as not measureable for oblate spheroids (aspect ratio < 1)
art_egg_database_filtered <- art_egg_database_raw %>% 
            filter(!(real_ar %in% c("0.5","1") & real_curv == "30")) %>% 
            mutate(im_curvature_deg = ifelse(real_ar %in% c("0.5","1"),NA,im_curvature_deg))

# Set up the reclicking dataframe
art_egg_database <- art_egg_database_filtered %>% rename(
            genus = cg,
            species = cs,
            image = im_status,
            # get measurements
            X1 = im_length_straight,
            X1px = im_length_straight_px,
            X2 = im_width,
            X2px = im_width_px,
            asym = im_asym,
            curv = im_curvature_deg
            # set up the final values
            ) %>% mutate(
            # calculating volume and aspect ratio from image measurements
            vol = ((4/3) * pi * (X1/2) * (X2/2)^2),
            ar = ifelse(!is.na(X1),
                        X1 / X2,
                        X1px / X2px),
            # calculate ellipticity as defined by Baker, Stoddard
            el = ar - 1,
            # calculate asymmetry as ratio - 1,
            asym = asym - 1,
            real_asym = real_asym - 1,
            # calculate differences
            diffar = real_ar - ar,
            diffcurv = real_curv - curv,
            diffasym = real_asym - asym,
            fracdiffar = diffar / real_ar,
            fracdiffcurv = diffcurv / real_curv,
            fracdiffasym = diffasym / real_asym,
            # log10 transformations
            logX1 = log10(X1),
            logX2 = log10(X2),
            logvol = log10(vol),
            logar = log10(ar),
            sqasym = sqrt(asym),
            sqcurv = sqrt(curv),
            logrealar = log10(real_ar),
            logrealasym = sqrt(real_asym),
            logrealcurv = sqrt(real_curv),
            real_ar = as.character(real_ar),
            real_curv = factor(real_curv,levels=c("0","30","120")),
            real_asym = as.character(real_asym)
            # Use these values to remove curvature for aspect ratios equal to or below 1
            ) %>% mutate(curv = ifelse(real_ar %in% c("0.5","1"),NA,curv)
            ) %>% mutate(sqcurv = ifelse(real_ar %in% c("0.5","1"),NA,sqcurv)
            # these final transformations select the main measurements used downstream
            ) %>% select(
                genus,species,
                X1,X2,logX1,logX2,ar,
                el,curv,asym,vol,logar,logvol,
                diffar,diffasym,diffcurv,
                fracdiffar,fracdiffasym,fracdiffcurv,
                real_ar,real_curv,real_asym,
                sqasym,sqcurv,logrealar,logrealasym,logrealcurv
            )

cols <- c("black","#696969","#C0C0C0","white")
shaps <- c(21,22,23,24)

gcurv <- ggplot(art_egg_database %>% filter(!(real_ar %in% c(0.5,1))),aes(x = diffcurv, y = factor(real_curv,levels=c("0","30","120")), shape = real_asym, fill = real_ar)) + 
        geom_jitter(size=6,alpha=0.9,color = "black") + 
        scale_shape_manual(values = shaps) + 
        scale_fill_manual(values = cols[c(1,4)]) + 
        xlab("curvature difference") + 
        ylab("curvature") + 
        guides(fill = guide_legend(title = "aspect ratio",override.aes = list(shape = 21)),
                shape = guide_legend(title = "asymmetry"))


gasym <- ggplot(art_egg_database,aes(x = diffasym, y = real_asym, shape = real_curv, fill = real_ar)) + 
        geom_point(size=6,alpha=0.9,color = "black") + 
        scale_shape_manual(values = shaps) + 
        scale_fill_manual(values = cols) + 
        xlab("asymmetry difference") + 
        ylab("asymmetry") + 
        guides(fill = guide_legend(title = "aspect ratio",override.aes = list(shape = 21)),
                shape = guide_legend(title = "curvature"))

gar<- ggplot(art_egg_database,aes(x = diffar, y = real_ar, fill = real_asym, shape = real_curv)) + 
        geom_point(size=6,alpha=0.9,color = "black")  + 
        scale_shape_manual(values = shaps) + 
        scale_fill_manual(values = cols) + 
        xlab("aspect ratio difference") + 
        ylab("aspect ratio") + 
        guides(fill = guide_legend(title = "asymmetry",override.aes = list(shape = 21)),
                shape = guide_legend(title = "curvature"))


pdf(file="Figure_art_reclicking.pdf",height=14,width=8,useDingbats=F)
print(grid.arrange(gcurv,gasym,gar,ncol=1))
dev.off()

table_art <- art_egg_database %>% group_by(real_ar,real_asym,real_curv) %>% 
    summarise(mean_diffar = mean(diffar,na.rm=T),
        mean_diffasym = mean(diffasym,na.rm=T),
        mean_diffcurv = mean(diffcurv,na.rm=T)) %>% 
    arrange(real_ar,real_asym,real_curv) %>% 
    select(real_ar,real_asym,real_curv,mean_diffar,mean_diffasym,mean_diffcurv) %>% 
    rename("aspect ratio" = real_ar,
        "asymmetry" = real_asym,
        "curvature" = real_curv,
        "mean aspect ratio discrepancy" = mean_diffar,
        "mean asymmetry discrepancy" = mean_diffasym,
        "mean curvature discrepancy" = mean_diffcurv)


print(xtable(table_art,digits = c(2,2,2,2,2,2,2),align = c("|l|","|l|","l|","l|","p{0.1\\textwidth}|","p{0.1\\textwidth}|","p{0.1\\textwidth}|")),include.rownames=F,file="art_table_latex.txt")


