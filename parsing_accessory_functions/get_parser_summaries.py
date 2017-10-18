### get_parser_summaries.py
# This code was developed by SH Church starting in Oct, 2016
# The purpose is to summarize the progress in parsing publications.

# This code is not ready for out of the box use for purposes outside
# of the creation of the insect egg database.

import sys
from os import listdir
import os.path
import re
from time import gmtime, strftime
import pandas

def get_unique_journals(current_dir):
	# PURPOSE: Opens a directory of .bib files and then 
	# counts the number of unique journal names in full set

	# ARGUMENT:
	# current_dir : string - The relative path to a directory

	# RETURNS: nothing

	all_journals = []
	for f in listdir(current_dir):
		fPath = os.path.join(current_dir,f)
		# Checks if path is a file and checks for pdf extension
		if fPath[-4:]=='.bib':
			print fPath
			with open(fPath, 'rb') as f:
				for line in f:
					tmp = re.match(r'journal = {(.*)}', line)
					if tmp:
						journal = tmp.group(1)
						if journal not in all_journals:
							all_journals.append(journal)
	print len(all_journals)

def get_ref_count(bib_path):
	# PURPOSE: Opens a .bib file and returns the number of refs

	# ARGUMENT:
	# bib_path : string - The relative path to a .bib file

	# RETURNS:
	# refs : integer - The number of references in the .bib file

	refs = 0
	with open(bib_path, 'rb') as f:
		for line in f:
			tmp = re.match(r'^@.*{', line)
			if tmp:
				refs += 1
	return refs

def get_done_bib_info(text_path):
	# PURPOSE: Calculates stats from a done_bid_ids.txt file
	# Prints a few interesting quantities

	# ARGUMENT:
	# text_path : string - The relative path to a results file
	# e.g. ".../done_bib_ids.txt"

	# RETURNS:
	# stats : list - A list of integers, for each quantity
	# that is described below.

	# If the path is not real, returns a list of dashes. 

	if not os.path.isfile(text_path):
		return ['-','-','-','-']

	# Examined refs: the number that we've tried to download and/or open
	ex_refs = 0

	# Not online refs: self-explanatory
	no_refs = 0

	# Keep refs: we have the PDF for these, but they have no egg dimension data
	kp_refs = 0

	# Good refs: we've opened the PDF, and they contained good data
	gd_refs = 0

	with open(text_path, 'rb') as f:
		for line in f:
			ex_refs+=1
			if "#KP" in line or '#kp' in line:
				kp_refs+=1
			if "#NO" in line or "#no" in line:
				no_refs+=1
	gd_refs = ex_refs - no_refs - kp_refs

	return [ex_refs, no_refs, kp_refs, gd_refs]

def get_parsed_info(text_path):
	# PURPOSE: Calculates stats from a egg dict file
	# Prints the number of unique genera and species

	# ARGUMENT:
	# text_path : string - The relative path to a results file
	# e.g. ".../Hemiptera_all.txt"

	# RETURNS:
	# stats : list - A list of integers, for each quantity
	# that is described below.

	# If the path is not real, returns a list of dashes.

	if not os.path.isfile(text_path):
		return ['-','-','-']

	# Total number of entries
	entries = 0

	# Unique genus strings
	unique_genera = []

	# Unique genus+species pairs
	unique_species = []

	with open(text_path, 'rb') as f:
		for line in f:
			entries+=1
			# turns each line back into a dict
			# will throw an error here if the line is not a dict
			try:
				d = eval(line)
			except SyntaxError:
				print '\nWARNING! Line will not neatly turn into a dict.\n'
				print 'Problem file:', text_path, '\n'
				print 'Problem line:', line
			# check if dict has a 'g' key
			if 'g' in d:
				# if the genus isn't in genera list, append it
				if str(d['g']) not in unique_genera:
					unique_genera.append(str(d['g']))
				# if the genus+species string isn't in species
				# list, append it to the list
				if 's' in d:
					tmp = str(d['g'])+str(d['s'])
					if tmp not in unique_species:
						unique_species.append(tmp)

	return [entries, len(unique_genera), len(unique_species)]

if __name__=='__main__':

	# Makes a list of the paths to all the .bib files in "bibs_by_order"
	bib_paths = []
	for f in listdir("bibs_by_order"):
		fPath = os.path.join("bibs_by_order",f)
		# Checks if path is a file and checks for bib extension
		if os.path.isfile(fPath) and fPath[-4:].lower() == '.bib':
			bib_paths.append(fPath)

	all_data = [['ORDER', 'TO_REFS', 'EX_REFS', 'NO_REFS', 'KP_REFS', 'GD_REFS', 'ENTRIES', 'U_GEN', 'U_SP']]

	total_entries = 0
	total_genera = 0
	total_sp = 0

	# Loop over all orders in "bibs_by_order"
	for bib in bib_paths:
		one_order_data = []
		tmp = re.match(r'bibs_by_order/(.*)\.bib', bib)

		# Check to make sure order path name fits the pattern
		if tmp:
			order = tmp.group(1)

			# Get the ref count
			one_order_data = [order, get_ref_count(bib)]

			# Get ref info from "done_bib_ids.txt"
			single_done_bib_path = 'data_by_order/' + order + '_data/done_bib_ids.txt'
			one_order_data += get_done_bib_info(single_done_bib_path)

			# Get parsed entry data from "[order]_all.txt"
			single_order_parsed_path = 'data_by_order/' + order + '_data/' + order + '_all.txt'
			one_order_data += get_parsed_info(single_order_parsed_path)

			# Add to the running totals for entries and species
			if os.path.isfile(single_order_parsed_path):
				total_sp += one_order_data[-1]
				total_genera += one_order_data[-2]
				total_entries += one_order_data[-3]
		else:
			print "Regular expression not matched"
			break

		# Turn all entries into strings for easier display
		all_data.append(map(str,one_order_data))

	print ''
	print 'Parsing progress as of', strftime("%Y-%m-%d %H:%M:%S", gmtime())
	print ''
	
	# Prints the data in neatly formatted columns
	widths = [max(map(len, col)) for col in zip(*all_data)]
	for row in all_data:
		print "  ".join((val.ljust(width) for val, width in zip(row, widths)))

	print ''
	print 'Total entries:', total_entries
	print 'Total unique genera:', total_genera
	print 'Total unique species:', total_sp
	print ''

	dt = pandas.DataFrame(all_data)

	csv = pandas.DataFrame.to_csv(dt,sep="\t",header=False,index=False)   
	outfile_csv = "parser_summaries_report.csv"
	with open(outfile_csv,'w') as of_csv:
	   of_csv.write(csv)