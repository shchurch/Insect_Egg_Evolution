### make_freq_table.py
### This code was written by Seth Donoughe August, 2016.

# The purpose is to build a word frequency table which is used as a 
# reference for searching for scientific names in a PDF.

import functions_for_parser as ffp

if __name__=='__main__':

	# path to directory of PDFs
	currentPath = '/Users/seth/Desktop/papers'

	filesOnly = ffp.findPdfPaths(currentPath)

	# Make a new dictionary:
	wordList = ffp.makeDictionary(filesOnly)
	# Make a text file where we will write dictionary
	f = open('wordFreq-v6-813.txt', 'wb')
	f.write(str(wordList))