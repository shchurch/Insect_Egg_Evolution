### concatenate_new.py
### This code was written by SH Church starting December 2016
# The purpose of this code is to concatenate the new entries into 
# the insect egg database following additional parsing. The final version
# of this code used to create the 2017 insect egg database
# identified new entries as lines not present in a hardcoded file,
# egg_database_raw.txt, from a pool of candidates in the folder
# data_by_order/*_all.txt. New entries were assigned an entry ID,
# the order information was populated according to the file of origin,
# and the tentative new database was saved as tmp_egg_database.txt.

# This code is not ready to be used out of the box for any purposes other than
# the insect egg database.

import sys
import re

# Initialize the expected list of paths in data_by_order
order_list = ('Zoraptera','WO_JV','Trichoptera','Thysanoptera','Strepsiptera','Siphonaptera','Raphidioptera','Psocoptera','Plecoptera','Phthiraptera','Phasmatodea','Orthoptera','Odonata','NO_Aug9','Neuroptera','Megaloptera','Mecoptera','Mantophasmatodea','Lepidoptera','July26_2017','July19_2017','Isoptera','Hymenoptera','Hemiptera','Grylloblattodea','Ephemeroptera','Embioptera','Diptera','Dermaptera','Coleoptera','Canadian_Zootaxa','Blattodea','Aug28_2017')

# START ALL ORDERS FILE
# this will be commented out after the first time

### candidate_lines = []
### for order in order_list:
### 	candidate_file = "data_by_order/" + order + "_data/" + order + "_all.txt"
### 	with open(candidate_file) as candidate_file:
### 		for line in candidate_file.read().splitlines():
### 			match_line = re.match('^{.*}$',line)
### 			if match_line:
### 				update_line = re.sub('}',', \'order1\': \'' + order + '\'}',line)
### 				candidate_lines.append(update_line)
### 
### all_file = "egg_database_raw.txt"
### with open(all_file, "w") as all_file:
### 	for c in candidate_lines:
### 		all_file.write(c + "\n")

# UPDATE ALL ORDERS FILE
all_file = "egg_database_raw.txt"
 
previous_lines = []
with open(all_file) as all_file:
    previous_lines = all_file.read().splitlines() 


from collections import Counter
all_lines1 = previous_lines
[k for k,v in Counter(all_lines1).items() if v>1]


all_lines2 = []

new_lines = []
for order in order_list:
	candidate_file = "data_by_order/" + order + "_data/" + order + "_all.txt"
	with open(candidate_file) as candidate_file:
		for line in candidate_file.read().splitlines():
			match_line = re.match('^{.*}$',line)
			if match_line:
				update_line = re.sub('}',', \'order1\': \'' + order + '\'}',line)
				all_lines2.append(update_line)
				if update_line not in previous_lines:
					new_lines.append(update_line)

for i in xrange(len(all_lines2)):
	if all_lines2[i] not in all_lines2:
		print "here\n"
		print all_lines2[i]

both_file = "egg_database_raw.txt"
with open(both_file, "w") as both_file:
	for p in previous_lines:
		both_file.write(p + "\n")
	for l in new_lines:
		both_file.write(l + "\n")

used_IDs = []
eggdb = []
with open("egg_database.txt") as f:
	for line in f:
		ln = eval(line)
		used = ln['ID']
		used_IDs.append(used)
		eggdb.append(ln)

### these are previous problem databases which no longer exist
#with open("problems_taxonomy.txt") as f:
#	for line in f:
#		ln = eval(line)
#		used = ln['ID']
#		used_IDs.append(used)
#
#with open("problems_dbl_chk.txt") as f:
#	for line in f:
#		ln = eval(line)
#		used = ln['ID']
#		used_IDs.append(used)

IDlist = range(int(max(used_IDs)), 1000000)
i = 0

updated = []
for line in new_lines:
	ln = eval(line)
	while(True):
		if str(IDlist[i]) in used_IDs:
			i = i + 1
		else:
			new_id = str(IDlist[i])
			used_IDs.append(new_id)
			break
	ln['ID'] = new_id
	updated.append(ln)

with open("tmp_new_egg_database.txt",'w') as out:
	for k in xrange(len(eggdb)):
		print >>out, eggdb[k]
	for j in xrange(len(updated)):
		print >>out, updated[j]
		