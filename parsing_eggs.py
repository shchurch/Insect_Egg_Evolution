### parsing_eggs.py
# This code was developed by 
# Sam Church, Jordan Hoffman, and Seth Donoughe
# Starting April, 2016

# The purpose is to manage the creation of an egg
# size database, processing bibtex style references
# to a python database of dictionaries with accompanying directory of pdfs

# The version of this code used to create the insect egg database
# includes a completely hardcoded dictionary of egg-specific terms.
# It therefore would not be very useful as an out-of the box tool 
# for other types of data collection.
# However, the structure of a program - including the use of a custom dictionary
# of hot-keys, the automatic parser for scientific names, and the management
# of entries as python dictionaries - could be applied to a variety of tasks.

# Future versions of this program are planned which would be generalizable
# to many types of literature-based phenotyping projects. 
# Changes that will be made include reading in the custom dictionary from a 
# separate configuration table the user could edit, and removing all references
# to insect order and local insect egg database directories.

import pyperclip
import time
import sys, getopt
import termios
import contextlib
import numpy as np
import webbrowser
import re
import sys
import os.path
import string
import functions_for_parser as ffp
from time import gmtime, strftime
#import speech_recognition as sr

def get_args(argv):
    inputfile = '' #global reference file
    outputfile = '' #global output file for results from parsing
    pdfdir = '' #gloabl directory to PDF repository
    restart = '' #global restart flag
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'hb:o:p:r', ['help','bib_file=','out_file=','pdf_directory=','restart'])
    except getopt.GetoptError:
        #usage()
        sys.exit(2)
   
    for opt, arg in opts:
        #help is broken currently
        if opt in ('-h', '--help'):
            usage()
            sys.exit(2)
        elif opt in ('-b', '--bib_file'):
            inputfile = arg
        elif opt in ('-o', '--out_file'):
            outputfile = arg
        elif opt in ('-p', '--pdf_directory'):
            pdfdir = arg
        elif opt in ('-r', '--restart'):
            restart = 1
        else:
            #usage()
            sys.exit(2)
    return inputfile,outputfile,pdfdir,restart

#PURPOSE: write out information on this session of parsing
def write_run_info(info):
    with open(info, 'w') as f:
        f.write("\n        ..ooOO EGGS OOoo..        \n")
        f.write(strftime("\n  Start time : %Y-%m-%d %H:%M:%S", gmtime()))
        if(restart):
            f.write("\n * restarting from previous run")
            f.write("\n * reading in bibtex file:" + inputfile)
        if os.path.exists(pdfdir):
            f.write("\n * pdf directory: " + pdfdir)
        else:
            f.write("\n * no pdf library specied")

# PURPOSE: interpret key strokes as input (?)
# Jordan Hoffman is the pro for this function
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

# PURPOSE: print the parsing dictionary of hotkeys
# can be called by the user at any time
def print_dictionary():
    print '=== ------------------------ ==='
    print '===      ~ DICTIONARY ~      ==='
    print ' w = open google scholar, queue title in clipboard'
    print ' b = print bib id, queue bib id in clipboard'
    print ' g = genus'
    print ' s = species'
    print ' i = page of pdf with image'
    print ' p = point estimate mode'
    print '   l = length, point estiamte'
    print '   w = width, point estimate'
    print '   b = breadth, point estimate'
    print ' a = average mode'
    print '   l = length, average'
    print '   w = width, average'
    print '   b = breadth, average'
    print ' d = standard deviation mode'
    print '   l = length, standard deviation'
    print '   w = width, standard deviation'
    print '   b = breadth, standard deviation'
    print ' m = minimum mode'
    print '   l = length, minimum in range'
    print '   w = width, minimum in range'
    print '   b = breadth, minimum in range'
    print ' x = maximum mode'
    print '   l = length, maximum in range'
    print '   w = width, maximum in range'
    print '   b = breadth, maximum in range'
    print ' e = edit mode for manual entries'
    print ' # = skip current bibtex ID without entry'
    print ' $ = displays current entry'
    print ' % = deletes last entry - doesn\'t clear current'
    print ' ^ = go back one bibtex entry'
    print ' T = perform PDF mining tools'
    print ' q = quit - appends data and break'
    print ' n = next entry - appends data'
    print ' N = next Bibtex id - appends data'
    print ' P = print this dictionary'
    print '=== ------------------------ ==='

# PURPOSE: Allow for manual entry, rather than clipboard entry
# Each ch=sys.stdin.read(1) takes keystroke
# Then var = raw_input waits for typing
# Returns the keystroke used and entry
def edit():
    print 'Edit manually [your edit will not appear as you type]'
    ch = sys.stdin.read(1)
    ch2 = ()
    var = ()
    global egg_lib #pdf library
    global wordList #word frequency table
    if ch=='s':
        var = raw_input("Please enter new species: ")
    elif ch=='g':
        var = raw_input("Please enter new genus: ")
    elif ch=='i':
        var = raw_input("Please enter new image page number: ")
    elif ch=='p':
        print 'point mode'
        ch2,var = edit_len_wid()
    elif ch=='a':
        print 'average mode'
        ch2,var = edit_len_wid()
    elif ch=='d':
        print 'standard deviation mode'
        ch2,var = edit_len_wid()
    elif ch=='m':
        print 'minimum mode'
        ch2,var = edit_len_wid()
    elif ch=='x':
        print 'maximum mode'
        ch2,var = edit_len_wid()
    else:
        var = raw_input("Enter whatever you want: ")
        print '\nwacky entry is',var
    return ch,ch2,var

#PURPOSE: Process bibtex file, build entries from keystrokes
# Program starts, loads bib tex file, then waits for entry
# Each ffp.cleanUpString(pyperclip.paste()) pulls from the clipboard
# Each ch = sys.stdin.read(1) reads a keystroke
# retuns local_data, a local dictionary with entries
def run():
    global done_ids,comt_ids,egg_lib #global varaibles in main
    print '\nHere we go! Starting parsing...'
    ALL = [] #array with results from parsing
    local_data = {} #current entry
    #Iterate through all bib IDS in the reference file
    i = 0 
    #Retrieve the IDs and titles from the references
    bib_ids,bib_titles = import_bibs()
    print "bib ID:",bib_ids[i]
    print "bib title:", bib_titles[i]
    #Read in typed input from the terminal
    with raw_mode(sys.stdin):
        try:    
            #Continue reading input until a break
            while(True):
                #Check if this bib ID has been parsed before
                while bib_ids[i] in done_ids or bib_ids[i] in comt_ids:
                    print bib_ids[i],'is already in done_bib_ids.txt'
                    if bib_ids[i] in comt_ids:
                        print 'with comment',comt_ids[bib_ids[i]]
                    i += 1
                    print 'skipping to next id:',bib_ids[i]
                    #Stop the program when all bib IDs are read
                    if i == (len(bib_ids) -1):
                        print 'no more bib ids'
                        i -= 1
                        break           
                ### official start of parsing an entry:
                local_data['b'] = bib_ids[i] #store bib ID
                #Tap into mac notification system and print entry
                #This function is in functions_for_parser.py (ffp)
                ffp.makeNotificationText(local_data)
                #Read in text entry from the terminal
                ch = sys.stdin.read(1)
                #Read in from the clipboard, and remove any non-ASCII, non-whitespace
                a = ffp.cleanUpString(pyperclip.paste())
                ### the following if statements use the internal dictionary of the program
                ### this can be printed with P
                ### future iterations of the program should introduce a flexible dictionary
                if ch == 'g':
                    #If genus and species are both on clipboad,
                    # separate them and put in local_data
                    tmp = re.match(r'(\w+)\s+(\w+\s*\w*)', a)
                    if tmp:
                        genus = tmp.group(1)
                        species = tmp.group(2)
                        print 'genus:', genus
                        print 'species:', species
                        local_data['s'] = species
                    #Otherwise just store as genus
                    else:
                        genus = a
                        print 'genus:',genus
                    local_data['g'] = genus #store genus
                #Continue to check the clipboard throughout the entry
                # is this necessary ? 
                a = ffp.cleanUpString(ffp.cleanUpString(pyperclip.paste())) #why twice?
                if ch == 's':
                    species = a
                    print 'species:',species
                    local_data['s'] = species #store species
                a = ffp.cleanUpString(pyperclip.paste())
                if ch == 'i':
                    image = a
                    print 'image page number:',image
                    local_data['i'] = image #store image page number
                a = ffp.cleanUpString(pyperclip.paste())
                ### these commands add an extra layer of dictionary 
                ### they parse the measurement data, which comes in various types
                ### average, range, or point estimate
                ### as well as various dimensions
                ### length, width, and breadth
                if ch == 'p':
                    print 'point mode'
                    #Function to process measurement data
                    ID,var1,var2 = key_len_wid(ch) #ID = l,w, or b 
                    #If two numbers passed, store first as length, second as width
                    if var2:
                        local_data['pl'] = var1
                        local_data['pw'] = var2
                        print 'point length:', var1
                        print 'point width:', var2
                    #Otherwise, use measurement ID as key
                    else:
                        local_data[ID] = var1
                a = ffp.cleanUpString(pyperclip.paste())
                if ch == 'a':
                    print 'average mode'
                    #Function to process measurement data
                    ID,var1,var2 = key_len_wid(ch)
                    #If two numbers passed, store first as average, second as sd
                    if var2:
                        local_data['a'+ID[1]] = var1
                        local_data['d'+ID[1]] = var2
                        if ID[1] == 'l':
                            dim = 'length'
                        elif ID[1] == 'w':
                            dim = 'width'
                        else:
                            dim = 'unexpected key'
                        print 'average '+dim+':', var1
                        print 'deviation '+dim+':', var2 
                    else:
                        local_data[ID] = var1
                a = ffp.cleanUpString(pyperclip.paste())
                if ch == 'd':
                    print 'standard deviation mode'
                    #Function to process measurement data
                    ID,var1,var2 = key_len_wid(ch)
                    #If two numbers passed, store first as average, second as sd
                    if var2:
                        local_data['a'+ID[1]] = var1
                        local_data['d'+ID[1]] = var2
                        if ID[1] == 'l':
                            dim = 'length'
                        elif ID[1] == 'w':
                            dim = 'width'
                        else:
                            dim = 'unexpected key'
                        print 'average '+dim+':', var1
                        print 'deviation '+dim+':', var2 
                    else:
                        local_data[ID] = var1
                a = ffp.cleanUpString(pyperclip.paste())
                if ch == 'm':
                    print 'minimum mode'
                    #Function to process measurement data
                    ID,var1,var2 = key_len_wid(ch)
                    #If two numbers passed, store first as min, second as max
                    if var2:
                        local_data['m'+ID[1]] = var1
                        local_data['x'+ID[1]] = var2
                        if ID[1] == 'l':
                            dim = 'length'
                        elif ID[1] == 'w':
                            dim = 'width'
                        else:
                            dim = 'unexpected key'
                        print 'minimum '+dim+':', var1
                        print 'maximum '+dim+':', var2 
                    else:
                        local_data[ID] = var1
                a = ffp.cleanUpString(pyperclip.paste())
                if ch == 'x':
                    print 'maximum mode'
                    #Function to process measurement data
                    ID,var1,var2 = key_len_wid(ch)
                    #If two numbers passed, store first as min, second as max
                    if var2:
                        local_data['m'+ID[1]] = var1
                        local_data['x'+ID[1]] = var2
                        if ID[1] == 'l':
                            dim = 'length'
                        elif ID[1] == 'w':
                            dim = 'width'
                        else:
                            dim = 'unexpected key'
                        print 'minimum '+dim+':', var1
                        print 'maximum '+dim+':', var2 
                    else:
                        local_data[ID] = var1
                a = ffp.cleanUpString(pyperclip.paste())
                ### this section allows the user to manually type a value
                ### rather than relying on the clipboard
                if ch == 'e':
                    #The edit function receives typed user entries
                    ch_mode,ch_lw,new = edit()
                    #If the entry has a length and width, take type of measurement and value
                    if(ch_lw):
                        ID = ch_mode + ch_lw
                    #Otherwise just take value
                    else:
                        ID = ch_mode
                    #The key is pulled from the output of edit as well
                    local_data[ID] = new
                    #This list contains the known outputs, prints a meaningful line for the user
                    if ID == 'g':
                        print '\ngenus:',local_data['g']
                    if ID == 's':
                        print '\nspecies:',local_data['s']
                    if ID == 'i':
                        print '\nimage page number:',local_data['i']
                    if ID == 'pl':
                        print '\nlength, point:',local_data['pl']
                    if ID == 'pw':
                        print '\nwidth, point:',local_data['pw']
                    if ID == 'al':
                        print '\nlength, average:',local_data['al']
                    if ID == 'aw':
                        print '\nwidth, average:',local_data['aw']
                    if ID == 'dl':
                        print '\nlength, standard deviation:',local_data['dl']
                    if ID == 'dw':
                        print '\nwidth, standard deviation:',local_data['dw']
                    if ID == 'ml':
                        print '\nlength, minimum in range:',local_data['ml']
                    if ID == 'mw':
                        print '\nwidth, minimum in range:',local_data['mw']
                    if ID == 'xl':
                        print '\nlength, maximum in range:',local_data['xl']
                    if ID == 'xw':
                        print '\nwidth, maximum in range:',local_data['xw']
                ### these commands allow the user to interact with the parser
                #Print the reference title and ID with 'b'
                if ch == 'b':
                    #Place the bib ID on the clipboard for easy saving of PDFs
                    pyperclip.copy(ffp.cleanUpString(local_data['b']))
                    print 'current bib ID is',local_data['b']
                    print 'bib title:', bib_titles[i]
                #Print the parser dictionary with 'P'
                if ch == 'P':
                    print_dictionary()
                #Open the pdf, or if not found, open the web database for searching
                if ch == 'w':
                    pdfPath = egg_lib + "/" + local_data['b'] + ".pdf"
                    web_path = "file://" + pdfPath
                    if os.path.exists(pdfPath):
                        #By default, the web browser is used to view PDFs
                        webbrowser.open(web_path)
                    else:
                        path = local_data['b'] #is this necessary ? 
                        if(bib_titles):
                            #Place the reference title in the clipboard for easy searching
                            paste_title = bib_titles[i]
                            pyperclip.copy(ffp.cleanUpString(paste_title))
                        else:
                            #Or if there is no title, print a warning
                            print 'no bib titles loaded'
                        try:
                            #These provide webbrowser options to the user, uncomment as desired
                            webbrowser.open("https://scholar-google-com.ezp-prod1.hul.harvard.edu/")
                            #webbrowser.open("https://scholar.google.com/")
                            #webbrowser.open("https://nrs.harvard.edu/urn-3:hul.eresource:gscholar")
                            #webbrowser.open("https://scholar.google.com/schhp?hl=en&as_sdt=0,22&inst=5823668996110809182")
                        except (KeyboardInterrupt, EOFError):
                            print 'weird character in title' #should this be somewhere else?
                            pass
                #Move to the next entry in the same reference with lowercase 'n'
                if ch == 'n':
                    #The function entry_check verifies that the finished entry has reasonable values
                    ld = entry_check(local_data)
                    #Once the entry is cleaned, pass it to the data array
                    local_data = ld
                    #Start a new entry in the same reference
                    print '\n * Enter Next Species Information * \n' # could be moved down a few lines for clarity
                    print 'current bib ID is',bib_ids[i]
                    ALL.append(local_data) 
                    write_tmp_outfile(ALL)
                    local_data = {}
                #Move to the next reference with uppercase 'N"'
                if ch == 'N':
                    #The function entry_check verifies that the finished entry has reasonable values
                    ld = entry_check(local_data)
                    #Once the entry is cleaned, pass it to the data array
                    local_data = ld
                    ALL.append(local_data)
                    #The skip flag is used in writing to the file done_bib_ids.txt, skipped entries are not written
                    skip_flag = ()
                    write_dones(bib_ids[i],skip_flag)
                    write_tmp_outfile(ALL)
                    local_data = {}
                    #Iterate to the next reference
                    i += 1
                    #Alert the user if there are no more references
                    if i == len(bib_ids):
                        print 'no more bib ids'
                        i -= 1
                    print '\n *** Next bibtex Entry *** \n'
                    print 'current bib ID is',bib_ids[i]
                #Skip the current reference ID without writing to done_bib_ids.txt
                if ch=='#':
                    print " press # again to skip and ignore"
                    print " or press anything else to comment" 
                    #If you skip (#) twice, leave no comment and don't print the ID
                    ch2 = sys.stdin.read(1)
                    #Otherwise, print the ID as done with a comment
                    if not ch2=='#':
                        #The comment is entered by typing in the terminal
                        comt = raw_input("Type comment here:")
                        write_dones(bib_ids[i],comt)
                    local_data={}
                    #Iterate to the next reference
                    i += 1
                    #Alert the user if there are no more references
                    if i == len(bib_ids):
                        print 'no more bib ids'
                        i -= 1
                    print "\nbib ID:",bib_ids[i]
                #Delete the previous entry with '%'s
                if ch== '%':
                    ALL.pop(-1)
                #Move to the previous reference with '^'
                if ch== '^':
                    i -= 1
                #Print the current entry with '$'
                if ch== '$':
                    print 'current entry'
                    print local_data
                #Try the automatic parsing functions with 'T'
                if ch== 'T':
                    #The function pdf_mine attempts to find the genus, species, and egg measurements automatically
                    genus,species,pdfPath,candidate_string = pdf_mine(bib_ids[i])
                    #If the PDF was available, set up the entry with the auto-parsed information
                    if(pdfPath):
                        #Load the line with egg description is loaded to the clipboard for easy searching
                        pyperclip.copy(ffp.cleanUpString(candidate_string))
                        #Open the PDF with the webbrowser
                        web_path = "file://" + pdfPath
                        webbrowser.open(web_path)
                        #Store the identified genus and species
                        #Egg measurements are not automatically stored, as they are usually wrong
                        local_data['g'] = genus
                        local_data['s'] = species
                #Quit the parser with 'q'
                if ch== 'q':
                    #Write out all entries as a potential safety measure
                    ALL.append(local_data)
                    write_tmp_outfile(ALL)
                    local_data = {}
                    break
        #If the interactive parser fails, pass and try again
        except (KeyboardInterrupt, EOFError):
            pass
    return ALL #all entries (all local_data)

#PURPOSE: attempt to automatically parse egg data from a PDF
#Using function_for_parser functions, attempt to get
#genus, species, and egg measurements automatically from the text
def pdf_mine(bib_id):
    # path to directory of PDFs
    global egg_lib #pdf library
    global wordList #word frequencey table

    pdfPath = egg_lib + "/" + bib_id + ".pdf"

    candidateString = ()
    genus = ()
    species = ()
    if(ffp.isPdfPath(pdfPath)):
        # Tries to identify genus and species
        print strftime("%Y-%m-%d %H:%M:%S", gmtime()), "Parsing PDF:", pdfPath, '\n'
        pdfText = ffp.pdfToString(pdfPath)
        genus, species = ffp.findGenusAndSpecies(pdfText, wordList)
        candidateString, topScore = ffp.findSentence(pdfText)
        print "Genus:", genus
        print "Species:", species, '\n'
        print "Best string had score", topScore, '\n'
        print candidateString, '\n'
    
        # Tries to identify the egg length automatically:
        eggLengthIndex = ffp.findKeyWord(candidateString,1)
        if eggLengthIndex:
            print "Best guess for egg length:"
            print ffp.sendBestNumber(eggLengthIndex,candidateString), '\n'

        # Tries to identify the egg width automatically:
        eggWidthIndex = ffp.findKeyWord(candidateString,2)
        if eggWidthIndex:
            print "Best guess for egg width:"
            print ffp.sendBestNumber(eggWidthIndex,candidateString), '\n'
        print "MINED."
    return genus,species,pdfPath,candidateString

def is_number(lentry):
    try:
        float(lentry)
        return True
    except ValueError:
        return False

#PURPOSE: check that known dictionary values are reasonable
#clean up an entry before storing it
#checks if numbers are numbers, removes trailing punctuation and spaces
def entry_check(local_data):
    ld = local_data
    listo = ('al','dl','pl','ml','xl','aw','dw','pw','mw','xw','i','ab','db','pb','mb','xb','g','s')
    for l in listo:
        lentry = ld.get(l)
        if (lentry):
            lentry_lead = re.match(r'^\s+(.+)',lentry)
            lentry_trail = re.search(r'(.+)\s+$',lentry)
            if lentry_lead:
                print lentry,"<= entry has leading non-alphanumeric"
                #These commands are used for more interactive entry cleaning
                #should they be removed ? 
                #print "Do you want to change to '" + lentry_lead.group(1) + "'? [y or n]"
                #ch = sys.stdin.read(1)
                #if ch == "y":
                ld[l] = lentry_lead.group(1)
                #    print "new entry = ",ld[l]
                #if ch == "n":
                #    print "no change"
                #sys.stdin.seek(0)
            if lentry_trail:
                print lentry,"<= entry has trailing non-alphanumeric"
                ld[l] = lentry_trail.group(1)
                #print "Do you want to change to '" + lentry_trail.group(1) + "'? [y or n]"
                #ch2 = sys.stdin.read(2)
                #if ch2 == "y":
                ld[l] = lentry_trail.group(1)
                #    print "new entry = ",ld[l]
                #if ch2 == "n":
                #    print "no change"
                #sys.stdin.seek(0)
    listw = ('g','s')
    for l in listw:
        lentry = ld.get(l)
        if (lentry):
            lentry_lead = re.match(r'^\W+(.+)',lentry)
            lentry_trail = re.search(r'(.+)\W+$',lentry)
            if lentry_lead:
                print lentry,"<= entry has leading non-alphanumeric"
                #print "Do you want to change to '" + lentry_lead.group(1) + "'? [y or n]"
                #ch = sys.stdin.read(1)
                #if ch == "y":
                ld[l] = lentry_lead.group(1)
                #    print "new entry = ",ld[l]
                #if ch == "n":
                #    print "no change"
                #sys.stdin.seek(0)
            if lentry_trail:
                print lentry,"<= entry has trailing non-alphanumeric"
                ld[l] = lentry_trail.group(1)
                #print "Do you want to change to '" + lentry_trail.group(1) + "'? [y or n]"
                #ch2 = sys.stdin.read(2)
                #if ch2 == "y":
                ld[l] = lentry_trail.group(1)
                #    print "new entry = ",ld[l]
                #if ch2 == "n":
                #    print "no change"
                #sys.stdin.seek(0)
    listn = ('al','dl','pl','ml','xl','aw','dw','pw','mw','xw','i','ab','db','pb','mb','xb')
    for l in listn:
        lentry = ld.get(l)
        if (lentry):
            lentry_comma = re.search(r'^\d+\,\d+$',lentry)
            if lentry_comma:
                lentry = string.replace(lentry, ',', '.')
                ld[l] = lentry
                print "new entry = ",ld[l]
            res = is_number(lentry)
            if res is False:
                while (res is False):
                    print lentry,"is not a number"
                    lentry = raw_input("Please type new entry:")
                    res = is_number(lentry)
                    print "\n"
                ld[l] = lentry
                print "new entry = ",ld[l]
    return ld

#PURPOSE: write out done_bib_ids.txt
def write_dones(bib_id,skip_flag):
    file = "done_bib_ids.txt"
    with open(file, 'a') as file:
        file.write(bib_id)
        if(skip_flag):
            file.write(" : #" + skip_flag)
        file.write("\n")

#PURPOSE: write out entries as they are parsed
def write_tmp_outfile(res):
    with open('tmp_out_file.txt', 'w') as file:
        for item in res:
            file.write("{}\n".format(item))

# PURPOSE: Mini routine used by each 'mode' to manually edit length and width
# Similar to above, but gets length and width and processes mode
# Mode can be 'point mode', 'average mode', etc.
def edit_len_wid():
    ch2 = sys.stdin.read(1)
    if ch2=='l':
        var = raw_input("Please enter new length: ")
    elif ch2=='w':
        var = raw_input("Please enter new width: ")
    elif ch2=='b':
        var = raw_input("Please enter new breadth: ")    
    else:
        var = raw_input("Enter whatever you want: ")
        print '\nwacky entry is',var
    return ch2,var

# PURPOSE: Mini routine used by each 'mode' to log length and width
# Similar to above, but gets length or width from the clipboard
# Mode can be 'point mode', 'average mode', etc.
def key_len_wid(ch):
    a = ffp.cleanUpString(pyperclip.paste())
    ch2 = sys.stdin.read(1)
    a = ffp.cleanUpString(pyperclip.paste())
    ch_lw = ch2
    ch_mode = ch
    ID = ch_mode + ch_lw
    a,b = separate_two_numbers(a)
    if not b:
        if ch2 == 'l':
            length = a
            print 'length:',length
        elif ch2 == 'w':
            width = a
            print 'width:',width
        elif ch2 == 'b':
            breadth = a
            print 'breadth:',breadth
        else:
            print 'unexpected key'
    return ID,a,b

# PURPOSE: Determine if a string contains two numbers,
# separated by characters that are neither numbers nor commas 
# or periods. If so, return those as first and second.
# Otherwise, return one number and None.
def separate_two_numbers(text):
    #tmp = re.match('\d+[.,]?\d*[^\d.,]+\d+[.,]?\d*',text)
    #tmp = re.match(r'\D*(\d+[.,]?\d*)[^\d.]+(\d+[.,]?\d*)\D*',text)
    tmp = re.match(r'\D*(\d+[.,]?\d*)[^\d.,]+(\d+[.,]?\d*)\D*',text)
    if tmp:
        first = tmp.group(1)
        first = string.replace(first, ',', '.')
        second = tmp.group(2)
        second = string.replace(second, ',', '.')
        return first, second
    else:
        return text, None

# PURPOSE: Processes bib file for bib ID and paper titles
# bib_ids = id for each bib entry
# bib_titles = correpsonding paper title
# checks that these are the same length or dies
def import_bibs():
    bib_file = inputfile
    bib_ids = []
    bib_titles = []
    f = open(bib_file,'r')
    for line in f:
        match_bib = re.match('\@.*\{(.*)\,',line)
        if match_bib:
            bib_ids.append(match_bib.group(1))
        match_title = re.match('\s*title\s*\=\s*\{+(.*?)\}+',line)
        if match_title:
            bib_titles.append(match_title.group(1))
    if len(bib_ids) != len(bib_titles):
        print 'bib tex file missing id or title'
        sys.exit()
    return bib_ids,bib_titles



if __name__ =='__main__':
    # Prints a header to start parsing
    print "\n        ..ooOO EGGS OOoo..        \n"
    done_file = "done_bib_ids.txt" #file with bib IDs that have been parsed
    #done_file = "lep_table.txt" #deprecated done file
    inputfile,outputfile,pdfdir,restart = get_args(sys.argv[1:])
    egg_lib = ''
    done_ids = [] #global list of done IDs
    comt_ids = {} #global list of comments on done IDs (stored as : #[XX])

    if (restart):
        #If restarting, look for done_bib_ids.txt
        if os.path.exists(done_file):
            print " * restarting from previous run"
            #Read in all the done IDs
            with open(done_file, "r") as file:
                for line in file.readlines():
                    match_comt = re.match('(.+)\s\:\\s(.*)$',line) #comment if exists
                    match_done = re.match('(.+)$',line) #done ID
                    if match_comt:
                        comt_ids[match_comt.group(1)] = match_comt.group(2)
                    if match_done:
                        done_ids.append(match_done.group(1))
        #If no done file, then can't restart
        else:
            print "ERROR: no restart file detected"
            sys.exit(2)
    else:
        #If not restarting, look for done_bib_ids.txt
        if (os.path.exists(done_file)): 
            #Stop the program before done file gets overwritten
            print " * Cannot overwrite done_bib_ids.txt"
            print "    ERROR: please use --restart or use a clean directory"
            sys.exit(2)
    
    #Hardcoded path to a word frequency table used in searching for scientific names
    dictPath = 'wordFreq-v5-813papers.txt'

    #Read in the word frequency table
    if os.path.exists(dictPath):
        print " * reading in word frequency table:", dictPath
        with open(dictPath, 'rb') as f:
            wordList = eval(f.read()) #frequency table
    else:
        #Warn if table doesn't exist
        print "\n    ERROR: no such frequency table file:",dictPath

    #Open file of references
    if os.path.exists(inputfile):
        print " * reading in bibtex file:",inputfile
    else:
        #Fail if it doesn't exist
        print "\n    ERROR: no such bibtex file:",inputfile
        sys.exit(2)
    #Open path to PDF repository
    if os.path.exists(pdfdir):
        egg_lib = pdfdir #global path to PDF repository
    else:
        #Warn if it doesn't exist, but don't fail
        print " * no pdf library specied"
        print "    pdf mining function 'T' not gonna work"
    #Set up output file
    if(outputfile):
        print " * output file:",outputfile
        info_file = "run_info_" + outputfile
        #Write a file with the information on this parsing session
        write_run_info(info_file)
        #Run the parser
        res = run() #res = results from parsing
        #Print the results of the parser when the program is stopped
        print res
        #Write the results to the output file
        with open(outputfile, 'w') as file:
            for item in res:
                file.write("{}\n".format(item))
    else:
        #If no output file, warn
        print " * WARNING: results not saved to an outfile"
        #Run the parser
        res = run()
        print res
