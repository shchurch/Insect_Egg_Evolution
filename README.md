### EGG_PIPELINE.SH

##This document was written by SH Church

##It describes the pipeline used to assemble the Insect Egg Database

##The usage of each internal program is described below with
# the purpose and arguments explained

##The version of this document used to assemble the Aug 2017
# draft of the database was written on Aug 30, 2017


### COLLECTING REFERENCES

##References were pulled from internet databases following an
# explicit criteria. These were stored as bibtex entries.

##The program set_bib_ids.py is used to introduce new
# bibs to the database, verifying that IDs are not overlapping.
# It collects all existing bib IDs from a directory 
# (hardcoded as bibs_by_order), and stores them
# in a list - bib_list.txt
#It takes as an input a bibtex file of new refs, and outputs
# a new file with the bib IDs updated.

#This program was only used for entries added after Oct 2016.
#Prior to this Mendeley was used to check for duplicates - 
# it often failed, so some duplicate IDs do exist in the database.

#Usage:
python set_bib_ids.py Diptera_tmp.bib Diptera_updated.bib
# [input] = Diptera_tmp.bib
#  file with new refs
# [output] = Diptera_updated.bib
#  file with corrected refs, to be added to the bibs directory
# additional written files = bib_list.txt   


### PARSING TEXT REFERENCES

##The program parsing_eggs.py is used to efficiently,
# consistently, and reproducibly extract descriptions 
# of insect eggs from pdfs.
#It is a combination of manual and automatic inputs, and 
# uses with an internal dictionary of hotkeys.

#Usage:
python parsing_eggs.py --b bibs_by_order/Diptera.bib --o Diptera_Aug28.txt --p ~/Dropbox/Sam_and_Seth/Papers/Diptera_bibids/
# --b = path to bibtex file with references to be parsed
# --o = output file name
# --p = path to pdf repository of references
# --r = optional restart flag
# additional written files = tmp_out.txt
#  temporary file storing new database entries
# done_bib_ids.txt
#  file with record of bib IDs already parsed

##The program get_parser_summaries.py calculates statistics
# on the state of parsing entries. It reads in the number
# of bibtex references in a hardcoded directory (bibs_by_order)
# and the number of database entries (data_by_order),
# and calculates the number left to parse.
#It also prints a report used in subsequent statistical analyses

#Usage:
python get_parser_summaries.py
# written files = parser_summaries_report.csv

##The program concatenate_new.py pulls new entries from a 
# hardcoded directory (data_by_order) and combines them into
# a raw database file and a database file with IDs added.
# New entries are identfied as those which are not 
# present in the raw database.

#Usage
python concatenate_new.py
# written files = egg_database_raw.txt
#  the raw database file, combining all entries in data_by_order
# tmp_new_egg_database.txt
#  a temporary file which can be converted to egg_database.txt
#  after verifying that data has not been lost
# To convert:
mv tmp_new_egg_database.txt egg_database.txt


### STILL TO CLEAN UP

### make_taxon_count_table.py
### make_freq_table.py
### functions_for_parser.py

#### THESE COMMANDS INTERACT WITH egg_database.py
#### Best practice is to create a temporary egg_database.py and then overwrite
### get taxonomy
python get_taxonomy.py egg_database.txt
########## - takes list of permitted orders

### clean egg_database
python clean_database.py egg_database.txt egg_database.txt ~/Dropbox/Sam_and_Seth/Papers/
### get bib info
python get_bibinfo.py egg_database.txt egg_database.txt
########### - need to change order to directory

### double check the entries
python double_check_entries.py egg_database.txt egg_database.txt

### convert the database to other formats
python convert_database.py egg_database.txt egg_database.csv

### analyze the database
R CMD BATCH analyze_data/R_code_build_dataframe.R