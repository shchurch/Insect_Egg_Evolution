### clean_database.py
# This code was developed by SH Church starting in Oct, 2016
# The purpose is to clean up known parsing errors in the egg database automatically.
# These include non numericals in numerical columns, non alphabet characters in names, etc.

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

def rinse(unwashed_entry):
	### change all integer values to strings
	ld = unwashed_entry
	ld_nonums = dict()
	for k, v in ld.items():
		if isinstance(v, tuple):
			v = str(v)
		try:
			int(v)
			new_v = str(v)
		except (ValueError):
			new_v = v
		ld_nonums[k] = new_v

	### remove empties
	ld_noempties = {k: v for k, v in ld_nonums.items() if v !=''}

	### get rid of leading spaces
	ld_nolead_space = dict()
	for k, v in ld_noempties.items():
		space = re.match(r'^\s+(.+)$',v)
		if space:
			print v, "<= entry has leading space"
			ld_nolead_space[k] = space.group(1)
		else:
			ld_nolead_space[k] = v
	ld_notrail_space = dict()
	for k, v in ld_nolead_space.items():
		space = re.match(r'^(.+)\s+$',v)
		if space:
			print v, "<= entry has trailing space"
			ld_notrail_space[k] = space.group(1)
		else:
			ld_notrail_space[k] = v

	rinsed = ld_notrail_space
	ld_notrail_space['cln'] = 'rnsd'
	return ld_notrail_space


def check_components(dat):
	d = dat
	with raw_mode(sys.stdin):
		g_flag = 0
		b_flag = 0
		missing_flag = 0
		double_flag = 0
		# Check for missing measurement and try to fix
		missing_flag = 0
		for letter in all_letters:
			avg_keys = ['a' + letter,'d' + letter]
			max_keys = ['m' + letter,'x' + letter]
			d_keys = d.keys()
			avg_len = len(set(avg_keys) & set(d_keys))
			max_len = len(set(max_keys) & set(d_keys))
			if avg_len == 1:
				missing_flag = 1
			if max_len == 1:
				missing_flag = 1
		if missing_flag == 1:
			d,missing_flag = fix_int_cols(d)
		if missing_flag > 0:
			d,missing_flag = fix_missing_cols(d)
		# Check for two measurements
		double_flag = 0
		for letter in all_letters:
			all_keys = ['a' + letter,'x' + letter,'p' + letter]
			d_keys = d.keys()
			all_len = len(set(all_keys) & set(d_keys))
			if(all_len > 1):
				d,double_flag = fix_double_cols(d)
		# Is there a bib ID?
		if('b' not in d):
			d,b_flag = fix_b_cols(d)
		# Is there a genus?
		if('g' not in d):
			d,g_flag = fix_g_cols(d)
		# If any of the above problems
		flags = missing_flag + double_flag + b_flag + g_flag
		if flags > 0:
			d['cln'] = ''
		# does it have any measurements at all?
		elif any(i in d for i in all_measurements):
			d['cln'] = 'comp'
		else:
			# If not, is there at least an image?
			if('i' in d):
				d['cln'] = 'comp'
			else:
				print "\nentry : \n"
				for k, v in d.iteritems():
					print "	",k," : ",v
				print "This entry has no measurements at all!"
				print "You should enter some and come back to cleaning"
				while True:
					ch = sys.stdin.read(1)
					#ch = 'y'
					if ch == 'T':
						pdf_open(d)
					elif ch == 'n':
						fix_flag = 1
						break				
					d['cln'] = 'rnsd'
	return d

def check_cols(dat):
	d = dat
	with raw_mode(sys.stdin):
		while True:
			extra_flag = 0 
			for i in d:
				if i not in all_expected_keys:
					tax_match = re.match(r'^tax_.*',i)
					bib_match = re.match(r'^bib_.*',i)
					im_match = re.match(r'^im_.*',i)
					prob_match = re.match(r'^problem.*',i)
					dbl_match = re.match(r'^dbl.*',i)
					if not(tax_match) and not(bib_match) and not (im_match) and not (prob_match) and not (dbl_match):
						extra_flag = 1
						bad_key = i
			if extra_flag == 1:
				print "\nentry : \n"
				for k, v in d.iteritems():
					print "	",k," : ",v				
				print "This entry has a wacky column " +bad_key+":"+d[bad_key]
				print "Do you want to delete it? [y/n/T]"
				ch = sys.stdin.read(1)
				if ch == 'T':
					pdf_open(d)
				if ch == 'y':
					d.pop(bad_key,"None")
					print "\nnew entry : \n"
					for k, v in d.iteritems():
						print "	",k," : ",v
				if ch == 'n':
					break
			else:
				break
		if extra_flag == 1:
			d['cln'] = 'comp'
		else:
			d['cln'] = 'cols'
	return d

def fix_int_cols(d):
	i = 0
	fix_flag = 1
	j_first_num = ''
	for j in d:
		if(j not in all_expected_keys):
			tax_match = re.match(r'^tax_.*',j)
			bib_match = re.match(r'^bib_.*',j)
			num_match = re.match(r'.*(\d)$',j)
			if not(tax_match) and not(bib_match):
				try:
					is_number(d[j])
					if num_match:
						j_first_num = num_match.group(1)
					int_k = j
					i = i + 1
				except ValueError:
					pass
	if i == 1:
		missing_k,num_missing = get_missing_k(d)
		if num_missing == 1:
			print "\nentry : \n"
			for k, v in d.iteritems():
				print "	",k," : ",v
			print "This entry looks like an easy fix, right?"
			print "Change this "  +int_k+":"+d[int_k]+ " to this " +missing_k+":"+j_first_num+d[int_k]+ " ? [y/n/e]"
			print "Type [e] to manual edit"
			while True:
				ch = sys.stdin.read(1)
				#ch = 'y'
				if ch == 'T':
					pdf_open(d)
				if ch == 'y':
					d[missing_k] = j_first_num+d[int_k]
					del d[int_k]
					fix_flag = 0
					print "\nnew entry : \n"
					for k, v in d.iteritems():
						print "	",k," : ",v
					break
				if ch == 'e':
					d[missing_k] = raw_input("Type new entry:")
					del d[int_k]
					fix_flag = 0
					print "\nnew entry : \n"
					for k, v in d.iteritems():
						print "	",k," : ",v					
					break
				elif ch == 'n':
					fix_flag = 1
					break
		else:
			fix_flag = 1
	return d, fix_flag

def fix_missing_cols(d):
	missing_flag = 1
	for letter in all_letters:
		avg_keys = ['a' + letter,'d' + letter]
		max_keys = ['m' + letter,'x' + letter]
		d_keys = d.keys()
		avg_len = len(set(avg_keys) & set(d_keys))
		max_len = len(set(max_keys) & set(d_keys))
		if avg_len == 1:
			for i in avg_keys:
				if (i not in d_keys):
					print "\nentry : \n"
					for k, v in d.iteritems():
						print "	",k," : ",v
					print "this entry is missing : ", i
					print "You can open the bibtex file with [T]"
					print "Do you want to change the value? [y/n]"
					print "or press p to change to point measurement"
					while True:
						ch = sys.stdin.read(1)
						if ch == 'T':
							pdf_open(d)
						if ch == 'y':
							missing_flag = 0
							d[i] = raw_input("Type new entry: ")
							print "\nnew entry : \n"
							for k, v in d.iteritems():
								print "	",k," : ",v							
							break
						if ch == 'p':
							missing_flag = 0
							for j in avg_keys:
								if j != i:
									old_key = j
							d['p' + letter] = d[old_key]
							del d[old_key]
							print "\nnew entry : \n"
							for k, v in d.iteritems():
								print "	",k," : ",v							
							break
						if ch == 'n':
							missing_flag = 1
							break
		if max_len == 1:
			for i in max_keys:
				if (i not in d_keys):
					print "\nentry : \n"
					for k, v in d.iteritems():
						print "	",k," : ",v
					print "this entry is missing : ", i
					print "You can open the bibtex file with [T]"
					print "Do you want to change the value? [y/n]"
					print "or press p to change to point measurement"
					while True:
						ch = sys.stdin.read(1)
						if ch == 'T':
							pdf_open(d)
						if ch == 'y':
							missing_flag = 0
							d[i] = raw_input("Type new entry: ")
							print "\nnew entry : \n"
							for k, v in d.iteritems():
								print "	",k," : ",v							
							break
						if ch == 'p':
							missing_flag = 0
							for j in max_keys:
								if j != i:
									old_key = j
							d['p' + letter] = d[old_key]
							del d[old_key]
							print "\nnew entry : \n"
							for k, v in d.iteritems():
								print "	",k," : ",v
							break
						if ch == 'n':
							missing_flag = 1
							break
	return d, missing_flag

def fix_double_cols(d):
	i = 0
	double_flag = 1
	for letter in all_letters:
		all_keys = ['a' + letter,'x' + letter,'p' + letter]
		d_keys = d.keys()
		combo = set(all_keys) & set(d_keys)
		if len(combo) > 1:
			print "\nentry : \n"
			for k, v in d.iteritems():
				print "	",k," : ",v
			print "this entry has both of these : ", combo
			print "You can open the bibtex file with [T]"	
			print "which would you like to remove? [a,d/m,x/p/n]"
			print "n means remove none of them"
			while True:
				ch = sys.stdin.read(1)
				#ch = 'n'
				if ch == 'T':
					pdf_open(d)
				if ch == 'a' or ch == 'd':
					d.pop('a' + letter,"None")
					d.pop('d' + letter,"None")
					double_flag = 0
					break
				elif ch == 'm' or ch == 'x':
					d.pop('m' + letter,"None")
					d.pop('x' + letter,"None")
					double_flag = 0
					print "\nnew entry : \n"
					for k, v in d.iteritems():
						print "	",k," : ",v
					break
				elif ch == 'p':
					d.pop('p' + letter,"None")
					double_flag = 0
					break
				elif ch == 'n':
					double_flag = 1
					break
	return d, double_flag

def fix_b_cols(d):
	b_flag = 1
	print "\nentry : \n"
	for k, v in d.iteritems():
		print "	",k," : ",v
	print "This entry is missing a bibtex! Thats no good"
	print "Do you know where to get one? [y/n]"
	ch = sys.stdin.read(1)
	if ch == 'y':
		d['b'] = raw_input("Type new bib ID: ")
		print "\nnew entry : \n"
		for k, v in d.iteritems():
			print "	",k," : ",v
		b_flag = 0
	else:
		b_flag = 1
	return d, b_flag

def fix_g_cols(d):
	g_flag = 1
	print "\nentry : \n"
	for k, v in d.iteritems():
		print "	",k," : ",v
	print "This entry is missing a genus!"
	print "You can open the bibtex file with [T]"
	print "Do you know where to get one? [y/n]"
	while True:
		ch = sys.stdin.read(1)
		if ch == 'T':
			pdf_open(d)
		if ch == 'y':
			d['g'] = raw_input("Type new genus: ")
			print "\nnew entry : \n"
			for k, v in d.iteritems():
				print "	",k," : ",v
			g_flag = 0
			break
		if ch == 'n':
			g_flag = 1
			break
	return d, g_flag

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

def get_missing_k(d):
	m = 0
	missing_key = ''
	for letter in ('l','w','b'):
		avg_key = ('a' + letter,'d' + letter)
		max_key = ('m' + letter,'x' + letter)
		if any(i in max_key for i in d):
			for i in max_key:
				if (i not in d):
					missing_key = i
					m = m + 1
		if any(i in avg_key for i in d):
			for i in avg_key:
				if (i not in d):
					missing_key = i
					m = m + 1
	return missing_key,m

def check_measurements(dat,letter):
	badguys = []
	goodguys = []
	for d in dat:
		point_key = 'p' + letter
		average_key = ('a' + letter,'d' + letter)
		max_key = ('m' + letter,'x' + letter)
		all_key = ('a' + letter,'m' + letter,'p' + letter)
		if point_key in d:
			other_keys = ('a' + letter,'d' + letter,'m' + letter,'x' + letter)
			if any(i in d for i in other_keys):
				badguys.append(d)
			else:
				goodguys.append(d)
		elif any(i in average_key for i in d):
			other_keys = ('m' + letter,'x' + letter)
			if any(i in d for i in other_keys):
				badguys.append(d)
			elif any(i not in d for i in average_key):
				badguys.append(d)
			else:
				goodguys.append(d)
		elif any(i in max_key for i in d):
			if any(i not in d for i in max_key):
				badguys.append(d)
			else:
				goodguys.append(d)
		else:
			goodguys.append(d)
	return badguys, goodguys

def is_number(lentry):
	try:
		float(lentry)
		return True
	except ValueError:
		return False

def check_vals(dat):
	ld = dat
	val_flag = 0
	with raw_mode(sys.stdin):
		for l in all_word_keys:
			lentry = ld.get(l)
			if (lentry):
				lentry_lead = re.match(r'^\W+(.+)',lentry)
				lentry_trail = re.search(r'(.+)\W+$',lentry)
				if lentry_lead:
					print lentry,"<= entry has leading non-alphanumeric"
					ld[l] = lentry_lead.group(1)
				if lentry_trail:
					print lentry,"<= entry has trailing non-alphanumeric"
					ld[l] = lentry_trail.group(1)
					ld[l] = lentry_trail.group(1)
		for l in all_num_keys:
			lentry = ld.get(l)
			if (lentry):
				res = is_number(lentry)
				if res is False:
					while (res is False):
						print lentry,"is not a number"
						lentry = raw_input("Please type new entry:")
						res = is_number(lentry)
						print "\n"
					ld[l] = lentry
					print "\nnew entry : \n"
					for k, v in ld.iteritems():
						print "	",k," : ",v
		for l in all_measurements:
			measur = ld.get(l)
			if (measur):
				if float(measur) > 30:
					print "\nentry : \n"
					for k, v in ld.iteritems():
						print "	",k," : ",v					
					print "\nWhoa, there is no way that "+l+" is "+measur
					print "You can open the bibtex file with [T]"
					print "Would you like to enter a new value? [y/n]"
					while True:
						ch = sys.stdin.read(1)
						if ch == 'T':
							pdf_open(ld)
						if ch == 'y':
							ld[l] = raw_input("Type new entry: ")
							print "\nnew entry : \n"
							for k, v in ld.iteritems():
								print "	",k," : ",v
							break
						if ch == 'n':
							val_flag = 1
							break
		for letter in all_letters:
			avg = ld.get('a'+letter)
			stdev = ld.get('d'+letter)
			if avg and stdev:
				if (float(stdev) >= float(avg)):
					print "\nnew entry : \n"
					for k, v in ld.iteritems():
						print "	",k," : ",v
					print "stdev is larger than average : "+stdev+">"+avg
					print "Want to switch them? [y/n]"
					while True:
						ch = sys.stdin.read(1)
						if ch == 'T':
							pdf_open(ld)
						if ch == 'y':
							ld['a'+letter] = stdev
							ld['d'+letter] = avg
							print "\nnew entry : \n"
							for k, v in ld.iteritems():
								print "	",k," : ",v
							break
						if ch == 'n':
							val_flag = 1
							break	
			maximum = ld.get('x'+letter)
			minimum = ld.get('m'+letter)
			if maximum and minimum:
				if (float(minimum) >= float(maximum)):
					print "\nnew entry : \n"
					for k, v in ld.iteritems():
						print "	",k," : ",v
					print "minimum is larger than maximum : "+minimum+" > "+maximum
					print "Want to switch them? [y/n]"
					while True:
						ch = sys.stdin.read(1)
						if ch == 'T':
							pdf_open(ld)
						if ch == 'y':
							ld['x'+letter] = minimum
							ld['m'+letter] = maximum
							print "\nnew entry : \n"
							for k, v in ld.iteritems():
								print "	",k," : ",v
							break
						if ch == 'n':
							val_flag = 1
							break	
	if val_flag < 1:
		ld['cln'] = 'sqky'
	return(ld)

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

def write_out_file(washed):
	with open(allfile,'w') as file:
		for i in washed:
			file.write("{}\n".format(i))

if __name__=='__main__':
	file = args.input

	unwashed = []

	pdf_dir = args.pdf_dir

	all_letters = ('l','w','b')
	all_word_keys = ('s','g','b','cg','cs','order1','cln')
	all_measurements = ('av','dv','mv','xv','pv','al','dl','pl','ml','xl','aw','dw','pw','mw','xw','ab','db','pb','mb','xb')
	all_num_keys = all_measurements + ('i','ID')
	all_extra_keys = ('!','W','pd')
	all_expected_keys = all_num_keys + all_word_keys + all_extra_keys

	with open(file) as f:
		for line in f:
			ln = eval(line)
			if(ln.get('cln')):
				unwashed.append(ln)
			else:
				ln['cln'] = ''
				unwashed.append(ln)

	allfile = args.output
	cleanfile = "only_clean_entries.txt"

	with open(cleanfile,'w') as cf:
		washed = unwashed
		i = 0
		for uw in unwashed:
			flag = 0
			if uw['cln'] == 'sqky':
				flag = 1
			if uw['cln'] == '':
				uw = rinse(uw)
			if uw['cln'] == 'rnsd':
				uw = check_components(uw)
			if uw['cln'] == 'comp':
				uw = check_cols(uw)
			if uw['cln'] == 'cols':
				uw = check_vals(uw)
			if uw['cln'] == 'sqky':
				washed[i] = uw
				#print >>cf,uw
			else:
				uw = rinse(uw)
				washed[i]
			#if flag == 0:
			#	write_out_file(washed)
			i = i + 1
		write_out_file(washed)







