#!/usr/bin/env python

### Created by Bruno de Medeiros (souzademedeiros@fas.harvard.edu), starting on aug-2016
### The purpose of the script is to count the number of species that we have for each named taxon
###    and obtain the number of decribed species from open tree of life
### This script uses takes as input a text file containing a list of dictionaries
### These dictionaries must have keys with taxonomic ranks in the form tax_rank.
### Example: tax_family, tax_suborder, etc.
### Results are written to a csv file, with each record containing:
### taxon, rank, parent taxon, # spp in our dataset, # spp in OTL


import argparse, requests, pandas as pd, sys
parser = argparse.ArgumentParser()
parser.add_argument('-i','--input', help = 'path to input file, containing one dictionary per line', default = './egg_database.txt')
args = parser.parse_args()



#This function uses Open Tree of Life API version 3(https://github.com/OpenTreeOfLife/germinator/wiki/Taxonomy-API-v3)
#Given a taxon name, it returns the number of tips in OTL, the corresponding node id and the taxon rank
def get_info_OTL(taxon):
    sys.stdout.write('Retrieving info for taxon ' + taxon + '\n')
    sys.stdout.flush()
    #first, look up Taxon name in OPL taxonomy to get
    r = requests.post('https://api.opentreeoflife.org/v3/tnrs/match_names',
                      data = {'names':[taxon, taxon], #needs to have at least two elements, otherwise error
                              'do_approximate_matching':False,
                              'context_name':'Arthropods'})
    ott_id = r.json()['results'][0]['matches'][0]['taxon']['ott_id']

    #now, retrieve number of species according to taxonomy alone
    r = requests.post('https://api.opentreeoflife.org/v3/taxonomy/taxon_info',
                    data = {'ott_id':long(ott_id),
                            'include_terminal_descendants':True,
                            'include_lineage':True})
    num_tips_tax = len(r.json()['terminal_descendants'])
    OTL_rank = r.json()['rank']
    for i in xrange(len(r.json()['lineage'])):
        parent = r.json()['lineage'][i]['name']
        if parent: #if name is empty, get the next one
            break


    #now, retrive node info in the consensus phylogeny
    r = requests.post('https://api.opentreeoflife.org/v3/tree_of_life/node_info',
                      data = {'ott_id':long(ott_id),
                              'include_lineage':True})
    try:
        OTL_node_id = r.json()['node_id']
        num_tips_tree = r.json()['num_tips']
    except:
        OTL_node_id = None
        num_tips_tree = None

    return {'taxon': taxon, 'OTL_rank':OTL_rank,'parent':parent, 'OTL_node_id':OTL_node_id, 'Ntips_phylo':num_tips_tree, 'Ntips_tax':num_tips_tax}

def get_num_in_DB(DB,taxon, rank):
    #DB is a list of dictionaries, taxon is the taxon name to be searched, rank is the taxon rank
    tax_key = 'tax_' + rank

    counter = 0

    for record in DB:
        try:
            if record[tax_key] == taxon:
                counter += 1
        except:
            pass

    return counter

if __name__ == "__main__":
    #read list of dictionaries
    with open(args.input, 'r') as infile:
        records = [eval(line) for line in infile]

    #loop through records, save all named taxa
    taxa = set()
    for i in xrange(len(records)):
        taxa.update([y for x,y in records[i].iteritems() if x.startswith('tax_') and
                                                            x in ['tax_family',
                                                                'tax_order',
                                                                'tax_genus',
                                                                'tax_class']
                                                            ])
                                                            ## these ones below are optional
                                                            #x != 'tax_matched' and
                                                            #x != 'tax_score' and
                                                            #x != 'tax_source' and
                                                            #x != 'tax_genus' and
                                                            #x != 'tax_species' and 
                                                            #x != 'tax_cg_ott_id' and 
                                                            #x != 'tax_cs_ott_id' and
                                                            #x != 'tax_cs_ott_id' and
                                                            #x != 'tax_cs_ott_id' and
                                                            #x != 'tax_cg_ncbi_id' and
                                                            #x != 'tax_cs_ncbi_id' and
                                                            #x != 'tax_infraclass' and 
                                                            #x != 'tax_infraorder' and 
                                                            #x != 'tax_parvorder' and 
                                                            #x != 'tax_subfamily' and 
                                                            #x != 'tax_subgenus' and 
                                                            #x != 'tax_suborder' and 
                                                            #x != 'tax_subtribe' and 
                                                            #x != 'tax_superfamily' and
                                                            #x != 'tax_superorder' and
                                                            #x != 'tax_tribe' and
                                                            #x != 'tax_species group' and 
                                                            #x != 'tax_species subgroup' and 
                                                            #x != 'tax_ott_version' and
                                                            #x != 'tax_matched_id_in_source' 
                                                            #x != 'tax_no rank' and #the lines below are temporary. records generated as of 26 aug 2016 do not have those keys
                                                            #x != 'tax_phylum' and
                                                            #x != 'tax_kingdom' and
                                                            #x != 'tax_domain' and
                                                            #x != 'tax_class' and
                                                            #x != 'tax_superclass' ])
    taxa = sorted(taxa)
    #loop through all taxa, make a list of dicts with OTL info
    taxa_DB = [get_info_OTL(taxon) for taxon in taxa]
    #now add our species counts
    for i in xrange(len(taxa_DB)):
        taxa_DB[i].update({'N_in_DB':get_num_in_DB(records,taxa_DB[i]['taxon'], taxa_DB[i]['OTL_rank'])})

    print taxa_DB

    #transform to a table and write
    taxa_table = pd.DataFrame(taxa_DB)

        #first, generate name of outfile
    with open('taxon_count.csv', 'w') as outfile:
        taxa_table.to_csv(outfile)
