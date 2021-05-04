This repository contains all the code required to reproduce the results of two manuscripts:

**A database of egg size and shape from more than 6,700 insect species.** Samuel H. Church, Seth Donoughe, Bruno A. S. de Medeiros, Cassandra G. Extavour
_Scientific Data_ 6, 104 (2019) <https://doi.org/10.1038/s41597-019-0049-y>; _bioRxiv_ <https://doi.org/10.1101/471953>

**Insect egg size and shape evolve with ecology but not developmental rate.** Samuel H. Church*, Seth Donoughe*, Bruno A. S. de Medeiros, Cassandra G. Extavour
_Nature_ 571, 58-62 (2019) <https://doi.org/10.1038/s41586-019-1302-4>; _bioRxiv_ <https://doi.org/10.1101/471946>

## Contents:

### analyze_data/

This directory contains all the R code and data tables required to reproduce the results in <https://doi.org/10.1101/471946>. 

### create_database/

This directory contains the code used in creating the database from published descriptions in the literature.

### phylogeny/

This directory contains the scripts used to generate the phylogeny for insect genera present in the egg database

### egg database, Nov 2018

`egg_database.tsv`

This tsv file is the version of the egg database used to analyze the data in November 2018, and is necessary for running the scripts in `analyze_data`. It contains all of the raw measurements extracted from publications. 

A version of the database containing calculated values (e.g. volume, aspect ratio) is included with the file `egg_database_final_values.tsv`.

The sources for the egg database are provided in `bibliography_egg_database.bib`

