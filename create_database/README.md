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

### set_bibids.py ###

Used to introduce new bibs to the database, verifying that IDs 
are not overlapping.
It collects all existing bib IDs from a directory 
(hardcoded as bibs_by_order), and stores them
in a list - bib_list.txt
It takes as an input a bibtex file of new refs, and outputs
a new file with the bib IDs updated.

This program was only used for entries added after Oct 2016.
Prior to this Mendeley was used to check for duplicates - 
it sometimes failed, so some duplicate IDs do exist in the database.

```
python set_bib_ids.py [input] [output]
[input] = file with new refs
[output] = file with corrected refs, to be added to the bibs directory
additional written files = bib_list.txt   
```

### PARSING TEXT REFERENCES ###

### parsing_eggs.py ###
 used to extract descriptions 
of insect eggs from pdfs.
It is a combination of manual and automatic inputs, and 
uses with an internal dictionary of hotkeys.

```
python parsing_eggs.py --b [input bibtex file] --o [output text file] --p [directory with pdfs]
--r = optional restart flag
additional written files = tmp_out.txt
 temporary file storing new database entries
done_bib_ids.txt
 file with record of bib IDs already parsed
```

### get_parser_summaries.py ###
Calculates statistics on the state of parsing entries. It reads in the number
of bibtex references in a hardcoded directory (bibs_by_order)
and the number of database entries (data_by_order),
and calculates the number left to parse.
It also prints a report used in subsequent statistical analyses

```
python get_parser_summaries.py
written files = parser_summaries_report.csv
```


### concatenate_new.py ###
Pulls new entries from a hardcoded directory (data_by_order) 
and combines them into a raw database file and a database 
file with IDs added. New entries are identfied as those which 
are not present in the raw database.

```
python concatenate_new.py
written files = egg_database_raw.txt
 the raw database file, combining all entries in data_by_order
tmp_new_egg_database.txt
 a temporary file which can be converted to egg_database.txt
 after verifying that data has not been lost
To convert:
mv tmp_new_egg_database.txt egg_database.txt
```

### STILL TO CLEAN UP ###

make_taxon_count_table.py
make_freq_table.py
functions_for_parser.py

clean egg_database
python clean_database.py egg_database.txt egg_database.txt ~/Dropbox/Sam_and_Seth/Papers/
get bib info
python get_bibinfo.py egg_database.txt egg_database.txt
- need to change order to directory

double check the entries
python double_check_entries.py egg_database.txt egg_database.txt

convert the database to other formats
python convert_database.py egg_database.txt egg_database.csv
