import pandas
import re
import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('input', help = 'path to input file, containing one dictionary per line')
parser.add_argument('output', help = 'path to output file, containing one dictionary per line')
args = parser.parse_args()


def is_number(lentry):
    try:
        float(lentry)
        return True
    except ValueError:
        return False

def entry_check(local_data):
    ld = local_data
    listo = ('al','dl','pl','ml','xl','aw','dw','pw','mw','xw','i','ab','db','pb','mb','xb','g','s')
    for l in listo:
        lentry = ld.get(l)
        if (lentry):
            try:
                lentry_lead = re.match(r'^\s+(.+)',lentry)
            except TypeError:
                print lentry
            lentry_trail = re.search(r'(.+)\s+$',lentry)
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

if __name__=='__main__':
    dat = []
    file = args.input
    with open(file) as f:
    	for line in f:
    		ln = entry_check(eval(line))
    		dat.append(ln)
    
    dt = pandas.DataFrame(dat)
    
    csv = pandas.DataFrame.to_csv(dt,sep="\t")   
    outfile_csv = args.output
    with open(outfile_csv,'w') as of_csv:
	   of_csv.write(csv)

    #jsn = pandas.DataFrame.to_json(dt)   
    #outfile_jsn = "egg_database_taxonomy.json"
    #with open(outfile_jsn,'w') as of_jsn:
        #of_jsn.write(jsn)

    #sql = pandas.DataFrame.to_sql(dt)   
    #outfile_sql = "egg_database_taxonomy.sql"
    #with open(outfile_sql,'w') as of_sql:
        #of_sql.write(sql)