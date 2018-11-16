# SOFTWARE

This document describes the software used to assemble the 
file egg_database.csv

The usage of each code file is described below with
the purpose and arguments explained

The version of this document used to assemble the Aug 2017
draft of the database was written on Aug 30, 2017

## COLLECTING REFERENCES

References were pulled from internet databases following an
explicit criteria (manuscript in preparation, 2018). These 
were stored as entries in a bibtex file.

The references which, following parsing, had egg measurements which are in the egg database are stored in the main directory of this project, under `bibliography_egg_database.bib`

The references which were collected and searched, but did not have egg measurements are in the file `examined_references_uncited.bib`.

### set_bibids.py ###

Used to introduce new bibs to the database, verifying that IDs 
are not overlapping. It collects all existing bib IDs from a directory 
(hardcoded as `bibs_by_order`), and stores them
in a list - `bib_list.txt`. It takes as an input a bibtex file of new refs, and outputs
a new file with the bib IDs updated.

This program was only used for entries added after Oct 2016.
Prior to this Mendeley was used to check for duplicates - 
it sometimes failed, so some duplicate IDs do exist in the database.

```
python set_bib_ids.py [input] [output]
```
`[input]` = file with new refs

`[output]` = file with corrected refs, to be added to the bibs directory

additional written files = `bib_list.txt`

## PARSING TEXT REFERENCES

### parsing_eggs.py ###
Used to extract descriptions of insect eggs from pdfs. It is a combination of manual and automatic inputs, and uses with an internal dictionary of hotkeys.

```
python parsing_eggs.py --b [input bibtex] --o [output] --p [directory with pdfs]
```
`--b` = reference file in bibtex format

`--o` = text file contaning list of python dictionaries with databse entries

`--p` = directory with pdfs, named according to bibtex IDs

`--r` = optional restart flag

additional written files = `tmp_out.txt`
temporary file storing new database entries

`done_bib_ids.txt`
file with record of bib IDs already parsed


### concatenate_new.py ###
Pulls new entries from a hardcoded directory (`data_by_order`) 
and combines them into a raw database file and a database 
file with IDs added. New entries are identfied as those which 
are not present in the raw database.

```
python concatenate_new.py
```
written files = `egg_database_raw.txt`
 the raw database file, combining all entries in `data_by_order`

`tmp_new_egg_database.txt`
 a temporary file which can be converted to egg_database.txt
 after verifying that data has not been lost

To replace database:
`mv tmp_new_egg_database.txt egg_database.txt`


### get_bib_info.py ###
Incorporates bibliographic information into the egg database by matching the dictionary key `b`
with the bib ID from a bibtex file, adds all fields as new keys in the dictionary.

`python get_bibinfo.py [input] [output]`

`input` = list of python dictionaries, with bib ID in key 'b'

`output` = list of python dictionaries, with bibliographic information in new keys

Bibliographic directory is hardcoded


### clean_egg_database.py ###
Automatically cleans known problems in the parsing output, such as trailing spaces, nonalphanumeric characters, etc.

`python clean_database.py [input] [output] [directory of pdfs]`

`input` = list of python dictionaries, with bib ID in key 'b'

`output` = list of python dictionaries, with bibliographic information in new keys

`pdf directory` = directory with pdfs, named according to bibtex IDs



### double_check_entries.py ###
Interactively allows the user to double check entries by ID by reopening the pdf
and presenting the collected data. Dictionary keys can be updated or marked as correct.

`python double_check.py [input] [output] [directory of pdfs]`
 
`input` = list of python dictionaries, to be cleaned

`output` = list of python dictionaries, updated

`pdf directory` = directory with pdfs, named according to bibtex IDs


### convert_database.py ###
Converts the format of the database. The output from `parsing_eggs.py` is a list of python dictionaries.
To analyze in R, this list is converted to a csv file. Other possible formats include JSON and sqlite.

`python convert_database.py [input] [output]`

`input` = text file with list of python dictionaries

`output` = csv file 


## ADDITIONAL SOFTWARE

### get_parser_summaries.py ###
Calculates statistics on the state of parsing entries. It reads in the number of bibtex references in a hardcoded directory (bibs_by_order) and the number of database entries (data_by_order), and calculates the number left to parse. It also prints a report used in subsequent statistical analyses

```
python get_parser_summaries.py
```

written files = parser_summaries_report.csv


### make_taxon_count_table.py ###
The purpose of the script is to count the number of species for each named taxon and obtain the number of described species from the Open Tree of LIfe. 

```
make_taxon_count_table.py -i [input]
```
`-i` = text file containing a list of dictionaries

### make_freq_table.py ###
Builds a word frequency table which is used as a reference for searching for scientific names in a PDF. Currently takes a hardcoded path to a directory of pdfs

```
make_freq_table.py
```

### functions_for_parser.py ###
The purpose of this code is to provide additional functions to the text parsing program, including automatic reading of text, searching 


