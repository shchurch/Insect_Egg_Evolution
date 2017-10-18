### set_bibids.py
# This code was developed by SH Church starting in Oct, 2016
# The purpose is to read a bibtex file and check that IDs are unique to the database
# If not, it corrects them
# It also prints a list of used IDs

# This code is not ready for out of the box use for purposes outside
# of the creation of the insect egg database.

import sys
import re
from os import listdir
import os.path
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('input', help = 'path to input directory, containing one dictionary per line')
parser.add_argument('output', help = 'path to output directory, containing one dictionary per line')
args = parser.parse_args()

#Read in list of used IDs
def get_bib_list():
	with open("bib_list.txt") as bl:
		bib_list = bl.read().splitlines()
	return bib_list #already used IDs

#Check if new ID is unique
def check_bib(bib_id):
	bib_list = get_bib_list() #already used IDs
	i = 0
	while(True):
		#Check if used already
		if bib_id in bib_list:
			#If so, add an extention and check again
			bib_id = bib_id + tag_list[i]
			#Iterate through extensions
			i = i+1
		else:
			#Add unique ID
			bibs_to_add.append(bib_id)
			break
	return bib_id

#Add new IDs to bib_list.txt
def update_bib_list(bib_id):
	outfile = "bib_list.txt"
	with open(outfile, 'a') as of:
		print >>of, bib_id

if __name__=='__main__':
	all_paths = [] #this stores all bibtex files
	all_list = [] #this stores all bibtex IDs

	#bibs_by_order is hardcoded bibtex storage location
	for f in listdir("bibs_by_order"):
		# Check if path is a file and checks for bib extension
		fPath = os.path.join("bibs_by_order",f)
		if os.path.isfile(fPath) and fPath[-4:].lower() == '.bib':
			#Add to path array
			all_paths.append(fPath)
	
	for bib in all_paths:
		with open(bib, 'rb') as f:
			for line in f:
				#Search for bib IDs (they start with @)
				bib_id = re.match(r'\@.*\{(.*)\,', line)
				if bib_id:
					#add to ID array
					all_list.append(bib_id.group(1))
	
	#Write out list of used IDs
	outfile = "bib_list.txt"
	with open(outfile, 'w') as of:
		for i in xrange(len(all_list)):
			print >>of, all_list[i]

	tag_list = ["a","b","c","d","e","f","g","a1","a2","a3","a4","a5"] #extensions to make unique IDs
	
	bib_file = args.input #input new reference file
	new_file = args.output #output corrected reference file
	bib_lines = [] #Array of all lines
	bib_ids = [] #Array of new IDs
	bib_list = [] #Array of corrected reference lines
	bibs_to_add = [] #Array of IDs added to bib_list.txt

	bib_article = "article" #Maybe unnecessary, just initializing variable
	with open(bib_file, 'rb') as f:
		for line in f:
			#Get new bibtex references
			bib_match = re.match(r'\@(.*)\{(.*)\,', line)
			if bib_match:
				bib_article = bib_match.group(1) #bib type
				bib_id = bib_match.group(2) #bib ID
				#Check if ID is unique
				bib_id = check_bib(bib_id) #new bib ID
				bib_line = "@" + bib_article + "{" + bib_id + ",\n" #corrected ID line
				#Add new ID to bib_list.txt
				update_bib_list(bib_id) 
				#Append new ID line
				bib_lines.append(bib_line)
			else:
				#Append any other line (not ID line)
				bib_lines.append(line)
		with open(new_file,'w') as nf:
			for i in xrange(len(bib_lines)):
				#Write new file
				nf.write(bib_lines[i])

	#Print to the terminal which bibs were added
	print "added bibs:"
	print bibs_to_add










