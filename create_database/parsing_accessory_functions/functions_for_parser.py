### functions_for_parser.py
### This code was created by Seth Donoughe September of 2016.

# The purpose of this code is to provide additional functions to the 
# text parsing program, including automatic reading of text, searching 
# for scientific names, and interacting with the Mac OS terminal.

# The subroutines in this program may be useful for other purposes
# outside of the insect egg database. However please note that
# the code currently works only with a Mac OS due to the notify
# subroutine.

from os import listdir
import os.path
from operator import itemgetter
import PyPDF2
import string
from time import gmtime, strftime
import sys
import signal
import unicodedata

# Used by the PDFMiner importer:
from cStringIO import StringIO
from pdfminer.pdfinterp import PDFResourceManager, PDFPageInterpreter
from pdfminer.converter import TextConverter
from pdfminer.layout import LAParams
from pdfminer.pdfpage import PDFPage

def handler(signum, frame):
	# Define a handler for the timeout
	# Written based on this:
	# http://stackoverflow.com/questions/492519/timeout-on-a-python-function-call
	print "PDF took too long to read"
	raise Exception("PDF reader timeout")

def findKeyWord(text, wordChoice):
	# Takes a chunk of text and returns the index of a word
	# of interest. If you want 'length' or a variant therein,
	# set wordChoice = 1. If you want 'width' set wordChoice = 2.
	
	if wordChoice == 1:
		keywords = ['length', 'long', 'Length', 'Long']
	elif wordChoice == 2:
		keywords = ['width', 'wide', 'Width', 'Wide', 'breadth', 'Breadth', 'diameter']
	else:
		print 'Set wordChoice to either 1 or 2.'
		return

	for word in keywords:
		if word in text:
			return text.index(word)

def sendBestNumber(wordIndex, text):
	# Takes an integer index and a string that contains
	# at least one digit.

	# wordIndex is the index of the last character in the
	# keyword (e.g. "width").

	# Walks forward and backward in the string from wordIndex.
	# Finds the closest number and returns it.

	if wordIndex > len(text)/2:
		maxIterations = wordIndex
	else:
		maxIterations = len(text)-wordIndex

	# i is the distance in the string from wordIndex
	i=1

	# alternately walk left and right away from the wordIndex
	# identify the closest number and store its index as numIndex

	while i < maxIterations:
		if wordIndex-i >= 0:
			if text[wordIndex-i].isdigit():
				numIndex = wordIndex-i
				minus = True
				break
		if wordIndex+i < len(text):
			if text[wordIndex+i].isdigit():
				numIndex = wordIndex+i
				minus = False
				break
		i+=1

	# number string to return
	numString = text[numIndex]

	# j is the distance in the string from
	# the beginning of the number
	j=1

	if minus:
		# walk left until a letter is found
		while True:
			if numIndex-j < 0:
				break
			else:
				if text[numIndex-j].isalpha():
					break
				numString=text[numIndex-j]+numString
				j+=1
	else:
		# walk right until a letter is found	
		while True:
			if numIndex+j >= len(text):
				break
			else:
				if text[numIndex+j].isalpha():
					break
				numString=numString+text[numIndex+j]
				j+=1

	numString = numString.strip()
	numString = numString.strip(string.punctuation)
	return numString.strip()

def makeDictionary(directory):
	# Takes a directory of PDFs, reads them in, separates words,
	# and then makes a dictionary with each word as a KEY, and its
	# frequency in the PDFs as its VALUE.

	freqDictionary = {}

	totalWords = 0
	for pdf in directory:
		# read in file, add to allText.
		print strftime("%Y-%m-%d %H:%M:%S", gmtime()), "Opening:", pdf
		text = pdfToString(pdf)
		print strftime("%Y-%m-%d %H:%M:%S", gmtime()), "PDF parsed into a string"
		for word in text.split():
			# makes lowercase and removes punctuation
			tmp = word.strip(string.punctuation).lower()
			totalWords += 1
			if tmp in freqDictionary:
				freqDictionary[tmp] += 1
			else:
				freqDictionary[tmp] = 1
	
	# Remove words that are only observed once
	# and convert counts to frequencies
	for key in freqDictionary.keys():
		if freqDictionary[key] == 1:
			del freqDictionary[key]
		else:
			freqDictionary[key]=float(freqDictionary[key])/totalWords

	# Print the total number of non-unique entries in the frequency table
	print "Total entries in frequency table:", len(freqDictionary.keys())
	return freqDictionary

def findGenusAndSpecies(text, freqDictionary):
	# Takes the first two pages of a PDF and generates a list of
	# candidate strings that could be the species name in the paper.

	speciesList = []
	latinEndings = ['a','e','i','m','n','o','r','s','x']

	# Makes a frequency dictionary of candidates
	singlePaperDictionary = {}
	total = 0
	for word in text.split():
		total += 1
		if len(word)>4 and word.isalpha() and word[0].islower():
			tmp = word.strip(string.punctuation).lower()
			if tmp in singlePaperDictionary:
				singlePaperDictionary[tmp] += 1
			else:
				singlePaperDictionary[tmp] = 1

	for key in singlePaperDictionary.keys():
		# checks that candidates occur >4 times and have a latin ending char
		if singlePaperDictionary[key] > 5 and key[-1] in latinEndings:
			currentFreq = float(singlePaperDictionary[key])/total
			# need to check whether this is a good way to check if a
			# key is in a dictionary:
			if key not in freqDictionary:
				relFreq = 100000
			else:
				generalFreq = freqDictionary[key]
				relFreq = currentFreq/generalFreq
			# print key, singlePaperDictionary[key], relFreq
			speciesList.append([key, singlePaperDictionary[key], relFreq])

	# sort by relative frequency, and then break ties with current frequency
	speciesList.sort(key=itemgetter(1))
	speciesList.sort(key=itemgetter(2))

	# bestSpeciesList is list of top ranked of species candidate strings
	bestSpeciesList = []
	for entry in speciesList[-1:-6:-1]:
		bestSpeciesList.append(entry[0])

	# returns '' and '' if there are no good candidate species
	if len(bestSpeciesList) == 0:
		return '', ''
	noPuncText = []
	# make a version of the text with punctuation stripped from word ends
	for word in text.split():
		noPuncText.append(word.strip(string.punctuation))

	# check to see if the word preceding a candidate species is capitalized
	# if so, that becomes the best guess for genus, and bestSpecies is updated
	genus = ''
	bestSpecies = bestSpeciesList[0]
	for species in bestSpeciesList:
		firstSpeciesIndex = noPuncText.index(species)
		if len(noPuncText[firstSpeciesIndex-1]) > 3:
			if noPuncText[firstSpeciesIndex-1][0].isupper():
				genus = noPuncText[firstSpeciesIndex-1]
				bestSpecies = species
				break

	return genus, bestSpecies

def findSentence(text):
	# Returns a candidate promising string

	# blockLength is length of string to consider at a time
	blockLength = 400

	# offSet is how much the sliding window shifts at a time
	offSet = 20

	# bestString is the current highest scoring string
	bestString = ''
	topScore = 0

	if len(text) < blockLength:
		return "Text too short to find sentence", 0

	for i in xrange(len(text)/offSet):
		sentence = text[i*offSet:i*offSet+blockLength]
		currentScore = scoreSentence(sentence)
		if currentScore > topScore:
			topScore = currentScore
			bestString = sentence

	return bestString, topScore

def scoreSentence(sentence):
	# Sentence is candidate string.
	# Returns an integer score based on how promising a sentence is.
	score = 0
	if 'egg' in sentence or 'Egg' in sentence:
		score += 5
	if 'mature oocyte' in sentence or 'Mature oocyte' in sentence:
		score += 2
	if 'long' in sentence or 'length' in sentence:
		score += 1
	if 'wide' in sentence or 'width' in sentence or 'breadth' in sentence:
		score += 1
	if 'curved' in sentence or 'chorion' in sentence:
		score += 1
	if 'mm' in sentence:
		score += 1
	if 'head' in sentence or 'thorax' in sentence:
		score -= 1

	# Check to see if there are any numbers in the string
	total = 0
	for ch in sentence:
		if str.isdigit(ch):
			total += 1
	if total > 3:
		score += 2

	return score

def isPdfPath(fPath):
	# Takes a path to a file as an argument
	# Returns True if the file is real and ends in '.pdf'
	# Returns False otherwise

	# Checks if path is a file and checks for pdf extension
	if os.path.isfile(fPath) and fPath[-4:].lower() == '.pdf':
		return True
	else:
		print "pdf path not real:",fPath
		return False

def findPdfPaths(currentPath):
	# Takes a path to a directory as an argument
	# Finds all the sub-directories and files in the directory
	# Checks that each file is real and ends in '.pdf'
	# Appends the paths for the pdfs to pdfPaths
	# Returns the list of paths

	pdfPaths = []
	for f in listdir(currentPath):
		fPath = os.path.join(currentPath,f)
	 	# Checks if path is a file and checks for pdf extension
		if isPdfPath(fPath):
			pdfPaths.append(fPath)
	return pdfPaths

def resave_PDF_subset_with_PyPDF2(pdf_path, page_list, out_path):
	# PURPOSE: Opens a .PDF file and a list of page numbers, then
	# re-saves it with a new name.

	# ARGUMENT:
	# pdf_path : string - Path to a .PDF file
	# page_list : list - List of integer page numbers to concatenate
	# out_path : string - Path to the output PDF file

	# RETURNS:
	# nothing

	# Register the signal function handler (for timeout of PDF processing)
	signal.signal(signal.SIGALRM, handler)

	# Define a 10 second timeout for PDF reading
	signal.alarm(20)

	try:
		# Makes a File object for the current pdf
		pdfFileObj = open(pdf_path, 'rb')

		# The PdfFileReader class takes a File object
		pdfReader = PyPDF2.PdfFileReader(pdfFileObj)

		# Check that the PDF reader can determine the number of pages
		try:
			pages = pdfReader.numPages
		except:
			return '##### Error parsing PDF #####'

		# Checks that the requested pages are in the PDF
		for i in page_list:
			if i not in range(pages):
				return '##### Page '+str(i)+' is not in PDF #####'

		# Initializes a PDF writer
		out_PDF = PyPDF2.PdfFileWriter()

		# Adds pages to the out_PDF
		for j in page_list:
			p = pdfReader.getPage(j)
			out_PDF.addPage(p)

		# Write to file
		with open(out_path, 'wb') as f:
		    out_PDF.write(f)
		    return 'success'

	except Exception, exc: 
		print exc
		text = "##### PDF took too long to read #####"

	# Cancel the timer if the PDF was parsed before timeout
	signal.alarm(0)

def getStringUsingPDFMiner(pdfPath, pageList=None):
	# Takes a path to a pdf file and a list of the desired page numbers.
	# Uses PDFMiner to open it and extract the text of the pages.
	# Returns a concatenated string.
	# If no page numbers are given, returns all pages.

	# This function was altered from:
	# https://www.binpress.com/tutorial/manipulating-pdfs-with-python/167

	# More details on its usage here:
	# http://www.unixuser.org/~euske/python/pdfminer/programming.html

	# Note: PDFMiner is a little smarter than PyPDF2 in one regard: if you
	# ask for a page that isn't there it doesn't throw an error

	# Register the signal function handler
	signal.signal(signal.SIGALRM, handler)

	# Define a 10 second timeout for PDF reading
	signal.alarm(10)

	try:
		if not pageList:
			pagenums = set()
		else:
			pagenums = set(pageList)

		output = StringIO()
		manager = PDFResourceManager()
		converter = TextConverter(manager, output, laparams=LAParams())
		interpreter = PDFPageInterpreter(manager, converter)

		infile = file(pdfPath, 'rb')
		try:
			for page in PDFPage.get_pages(infile, pagenums):
				interpreter.process_page(page)
			infile.close()
			converter.close()
			text = output.getvalue()
		except:
			text = '##### Error parsing PDF #####'
		output.close

	except Exception, exc: 
		print exc
		text = "##### PDF took too long to read #####"

	# Cancel the timer if the PDF was parsed before timeout
	signal.alarm(0)

	return text

def getStringUsingPyPDF2(pdfPath, pageList=None):
	# Takes a path to a pdf file and a list of the desired page numbers.
	# Uses PyPDF2 to open it and extract the text of the page
	# Returns a concatenated string.
	# If no page numbers are given, returns all pages.

	# Register the signal function handler
	signal.signal(signal.SIGALRM, handler)

	# Define a 10 second timeout for PDF reading
	signal.alarm(10)

	try:
		# Makes a File object for the current pdf
		pdfFileObj = open(pdfPath, 'rb')

		# The PdfFileReader class takes a File object
		pdfReader = PyPDF2.PdfFileReader(pdfFileObj)

		# The number of pages in the PDF
		try:
			pages = pdfReader.numPages
		except:
			return '##### Error parsing PDF #####'

		if not pageList:
			pageList = range(pages)

		# Initialize empty string to hold page contents
		text = ''

		# Extracts the text from the pages
		for page in pageList:	
			if page >= pages:
				text += "##### Requested page does not exist #####"
			else:
				pageObj = pdfReader.getPage(page)
				try:
					pageText = pageObj.extractText()
				except ValueError:
					pageText = "##### Error parsing page " + str(page) +' #####'
					print "There was an error parsing page", page
				text += pageText

	except Exception, exc: 
		print exc
		text = "##### PDF took too long to read #####"

	# Cancel the timer if the PDF was parsed before timeout
	signal.alarm(0)

	return text

def cleanUpString(text):
	# Takes a string and does three things:
	# 1) Tries to turn unicode special characters with ASCII equivalents
	# 2) If it can't, replaces unprintable characters with '?'.
	# 3) Replaces all types of whitespace with ' '.

	tmp = ''
	for ch in text:
		# replace unprintable or unicode characters
		if isinstance(ch, unicode):
			unicodedata.normalize('NFKD', ch).encode('ascii','replace')
		if ch not in string.printable:
			tmp += '?'
		# replace misc whitespace with a space
		elif ch in string.whitespace:
			tmp += ' '
		else:
			tmp += ch
	return str(tmp)

def shrinkWhiteSpace(text):
	# Takes a string and reduces any stretch of whitespace
	# to a single space.
	tmp = ''
	for i in xrange(len(text)):
		if text[i] not in string.whitespace:
			tmp += text[i]
		elif text[i] in string.whitespace \
		and text[i-1] not in string.whitespace:
			tmp += ' '
	return tmp

def removeOneSpace(text):
	# Takes a string that has an extra space inserted
	# between all the characters and removes one space
	# per stretch of whitespace.

	tmp = ''
	spaceRemovedAlready = False
	for i in xrange(len(text)):
		if text[i] not in string.whitespace:
			tmp += text[i]
			spaceRemovedAlready = False
		elif text[i] in string.whitespace and spaceRemovedAlready == True:
			tmp += text[i]
		elif text[i] in string.whitespace and spaceRemovedAlready == False:
			# Omit a space and flip toggle
			spaceRemovedAlready = True
	return tmp

def wordStats(text):
	# Takes a string, splits it into words
	# Returns mean word length and total word count

	listOfWords = text.split()
	totalLength = 0
	for i in listOfWords:
		totalLength += len(i)
	totalWords = len(listOfWords)
	if totalWords != 0:
		return float(totalLength)/totalWords, totalWords
	else:
		return 0, 0

def pdfToString(pdfPath):
	# Take a path to a PDF file, extracts text,
	# checks for quality, re-extracts if necessary,
	# cleans up the formatting and returns a
	# single string with all the cleaned up text.

	# use PDFMiner to get text for first two pages
	testText = getStringUsingPDFMiner(pdfPath, [0,1])
	mean, total = wordStats(testText)

	if mean < 2 and total > 150:
		# Remove the extra spaces that PDFMiner sometimes puts in
		finalRaw = getStringUsingPDFMiner(pdfPath)
		final = shrinkWhiteSpace(removeOneSpace(cleanUpString(finalRaw)))
	elif mean > 8 or total < 150:
		# Use PDF2 to get text instead
		finalRaw = getStringUsingPyPDF2(pdfPath)
		final = shrinkWhiteSpace(cleanUpString(finalRaw))
	else:
		# Use PDFMiner to get text
		finalRaw = getStringUsingPDFMiner(pdfPath)
		final = shrinkWhiteSpace(cleanUpString(finalRaw))

	mean, total = wordStats(final)
	if total < 150:
		print 'WARNING: Very little text'
	return final

def makeNotificationText(dictionary):
	# PURPOSE: take a dictionary and convert it into three strings that
	# can be sent to the notify function.

	# The title line always has the bibID, genus, and species on it,
	# while the other lines have the rest of the dictionary entries. 

	# This constant is the max number of characters that we would
	# like to be displayed on the second notification lines.
	# Any remaining info will fill (and run off) the edge of the third line.
	subtitle_length = 45

	title,subtitle,message = '','',''

	if 'b' in dictionary:
		title += 'b:'+str(dictionary['b'])
	if 'g' in dictionary:
		title += ', g:'+str(dictionary['g'])
	if 's' in dictionary:
		title += ', s:'+str(dictionary['s'])

	a = 'bsg'
	for key in dictionary.keys():
		if key not in a:
			tmp = str(key)+':'+str(dictionary[key])
			if len(subtitle+tmp) <= subtitle_length:
				if len(subtitle) == 0:
					subtitle += tmp
				else:
					subtitle += ', ' + tmp
			else:
				if len(message) == 0:
					message += tmp
				else:
					message += ', ' + tmp

	notify(title, subtitle, message)

def notify(title, subtitle, message):
	# PURPOSE: Send a notification to apple OS
	# Requires that "terminal-notifier" be installed first with Ruby

	# Code altered from here:
	# http://stackoverflow.com/questions/17651017/python-post-osx-notification

	# Formatting the string that tells the terminal notifier what to do
	t = '-title {!r}'.format(title)
	s = '-subtitle {!r}'.format(subtitle)
	m = '-message {!r}'.format(message)
	a = '{}[]'

	for char in a:
		t = t.replace(char,"")
		s = s.replace(char,"")
		m = m.replace(char,"")

	os.system('terminal-notifier {}'.format(' '.join([m, t, s])))

if __name__=='__main__':

	pass
