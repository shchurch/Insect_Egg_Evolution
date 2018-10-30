import sys
from os import listdir
import os.path
import re
import argparse
from time import gmtime, strftime

parser = argparse.ArgumentParser()
parser.add_argument('input', help = 'path to input file, containing one dictionary per line')
parser.add_argument('output', help = 'path to output file, containing one dictionary per line')
parser.add_argument('-o','--overwrite', help = 'use this flag to overwrite all taxonomic information', action = 'store_true')
args = parser.parse_args()

def build_bib_db(order,bib):
	bib_entry = dict()
	bib_list = []
	bib_file_entries = []
	bib_database = dict()
	with open(bib, 'rb') as f:
		for line in f:
			bib_id = re.match(r'\@.*\{(.*)\,', line)
			if bib_id:
				bib_list.append(bib_id.group(1))
			entry = re.match(r'(.*)\ \=\ \{(.*)}', line)
			if entry:
				dict1 = "bib_" + entry.group(1)
				dict2 = entry.group(2)
				dct = {dict1:dict2}
				bib_entry.update(dct)
			end_match = re.match(r'^\}$', line)
			if end_match:
				bib_file_entries.append(bib_entry)
				bib_entry = dict()
	if len(bib_file_entries) != len(bib_list):
		print "error parsing bib"
		sys.exit()
	for i in xrange(len(bib_list)):
		bib_dict = {bib_list[i]:bib_file_entries[i]}
		bib_database.update(bib_dict)
	return bib_database

if __name__=='__main__':

	# Makes a list of the paths to all the .bib files in "bibs_by_order"
	bib_paths = []
	for f in listdir("bibs_by_order"):
		fPath = os.path.join("bibs_by_order",f)
		# Checks if path is a file and checks for bib extension
		if os.path.isfile(fPath) and fPath[-4:].lower() == '.bib':
			bib_paths.append(fPath)

	# Loop over all orders in "bibs_by_order"
	bib_db = dict()
	for bib in bib_paths:
		one_order_data = []
		tmp = re.match(r'bibs_by_order/(.*)\.bib', bib)

		# Check to make sure order path name fits the pattern
		if tmp:
			order = tmp.group(1)

			bib_database = build_bib_db(order,bib)
			bib_db.update({order:bib_database})

	with open(args.input, 'r') as infile:
		records = [eval(line) for line in infile]

	new_records = []
   	with open(args.output,'w') as outfile:
   		for i in xrange(len(records)):
   			if 'bib_fullID' in records[i] and not args.overwrite:
   				print >>outfile, records[i]
   			else:
   				try:
   					order_in_question = records[i]['order1']
   					matched_order = bib_db[order_in_question]
					bib_in_question = records[i]['b']
   					matched_bib = matched_order[bib_in_question]
   					new_rec = records[i]
   					new_rec['bib_fullID'] = bib_in_question
   					new_rec.update(matched_bib)
   					new_records.append(new_rec)
   					print >>outfile, new_rec
   				except(KeyError):
   					sys.stdout.write("problem with bib_id : " + bib_in_question + " in order " + order_in_question + "\n")
   					continue
   			



