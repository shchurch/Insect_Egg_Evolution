### convert_database.py
# This code was developed by SH Church starting in Oct, 2016
# The purpose is to interactively double check entries in the database
# It currently accepts typed egg IDs or bibtex IDs, repons a pdfs
# and allows the user to update an entry.

# This code is not ready for out of the box use for purposes outside
# of the creation of the insect egg database.

import argparse
import re
import pandas
import os.path
import sys
import contextlib
import termios
import webbrowser
parser = argparse.ArgumentParser()
parser.add_argument('input', help = 'path to input file, containing one dictionary per line')
parser.add_argument('output', help = 'path to output file, containing one dictionary per line')
parser.add_argument('pdf_dir', help = 'path to pdf directory')
args = parser.parse_args()

def pdf_open(d):
	bib_id = d['b']
	order = d['order1']
	if bib_id and order:
		#pdf_directory = "/Users/samuelchurch/Dropbox/Sam_and_Seth/Papers/"
		pdf_directory = pdf_dir
		pdf_path = pdf_directory + order + "_bibids/" + bib_id + ".pdf"
		web_path = "file://" + pdf_path
		if os.path.exists(pdf_path):
			webbrowser.open(web_path)
		else:
			print "pdf path doesn't exist : ",pdf_path
	else:
		print "missing bib ID or order, can't open pdf"

def edit():
	ch = raw_input("What key would you like to edit? (or n):")
	if ch == "n":
		var = ""
	else: 
		var = raw_input("Enter new value: ")
	return ch,var

@contextlib.contextmanager
def raw_mode(file):
    old_attrs = termios.tcgetattr(file.fileno())
    new_attrs = old_attrs[:]
    new_attrs[3] = new_attrs[3] & ~(termios.ECHO | termios.ICANON)
    try:
        termios.tcsetattr(file.fileno(), termios.TCSADRAIN, new_attrs)
        yield
    finally:
        termios.tcsetattr(file.fileno(), termios.TCSADRAIN, old_attrs)

def chk_entry(r):
	pdf_open(r)
	with raw_mode(sys.stdin):
		while(True):
			print "Would you like to edit? [e,n,$]"
			ch = sys.stdin.read(1)
			ch2 = ()
			var = ()
			if ch == "$":
				print 'current entry'
				for k, v in r.iteritems():
					print "	",k," : ",v
			elif ch == "n":
				break
			elif ch == "e":
				ch2 = raw_input("\nWhich key?: ")
				var = raw_input("\nEnter new value: ")
				print '\nnew entry is',ch2," : ",var
				r[ch2] = var
	var = raw_input("Enter dbl_chk flag: ") 
	r["dbl_chk"] = var
	return r

def write_out_file(washed):
	with open(outfile,'w') as file:
		for i in washed:
			file.write("{}\n".format(i))

if __name__=='__main__':
	file = args.input

	raw = []

	pdf_dir = args.pdf_dir
	outfile = args.output

	with open(file) as f:
		for line in f:
			ln = eval(line)
			if(ln.get('cln')):
				raw.append(ln)
			else:
				ln['cln'] = ''
				raw.append(ln)

	pdf_dir = args.pdf_dir
	
	with open(outfile,'w') as cf:
		while(True):
			Bib_in_q = raw_input("Type bib or ID to search:")
			dblchk = raw
			i = 0
			for r in raw:
				flag = 0
				if r["ID"] == Bib_in_q:
					print "dbl_chk flag is : ",r.get('dbl_chk')
					r = chk_entry(r)
					dblchk[i] = r
				elif r["b"] == Bib_in_q:
					Taxon_in_q = raw_input("Type taxon to search:")
					tax = Taxon_in_q.replace("_"," ")
					if r.get("g"):
						name_g = r.get("g")
					else:
						name_g = ""
					if r.get("s"):
						name_s = r.get("s")
					else:
						name_s = ""
					name = name_g + " " + name_s
					if(r.get("tax_matched") == tax):
						print "dbl_chk flag is : ",r.get('dbl_chk')
						r = chk_entry(r)
						dblchk[i] = r
					elif(name == tax):
						print "dbl_chk flag is : ",r.get('dbl_chk')
						r = chk_entry(r)
						dblchk[i] = r
				else:
					flag = 1
				if flag == 0:
					write_out_file(dblchk)
				i = i + 1

