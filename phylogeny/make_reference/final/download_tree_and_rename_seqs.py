#!/usr/bin/env python

### Created by Bruno de Medeiros (souzademedeiros@fas.harvard.edu), starting on Nov-2016
### The purpose of the script is finalize files that will be used for phylogenetic placement
### Misof tree is downloaded from Open Tree of Life
### 18S and 28S sequences are renamed using the same standard for SILVA sequences:
### Accession|ott id for species|ott id for genus|ott id for family|ott version|family|subfamily|genus|taxon name
### genera not represented in sequences are pruned from the tree
### tree is output as newick, renamed sequences as fasta


from Bio import SeqIO
from pandas import DataFrame
import argparse, requests, dendropy


#first, read input and save information in a pandas dataframe
parser = argparse.ArgumentParser()
parser.add_argument('-1', '--g18S', help = 'path to input 18S fasta file', default = '18S_to_rename.fasta')
parser.add_argument('-2', '--g28S', help = 'path to input 18S fasta file', default = '28S_to_rename.fasta')
parser.add_argument('-t','--table',help = 'path to table translating acession numbers to species names', default = 'reference_seqs_info.csv')
args = parser.parse_args()
#args = parser.parse_args([])

#read sequence files and reference table
g18S = [record for record in SeqIO.parse(args.g18S,'fasta')]
g28S = [record for record in SeqIO.parse(args.g28S,'fasta')]
reference = DataFrame.from_csv(args.table)



#create a dictionary acession numbers to species names
accessions = []
for seq in g18S:
    try:
        accession = seq.id.split('|')[0].split(':')[1]
    except IndexError:
        accession = seq.id.split('_')[0]
    accessions.append(accession)
for seq in g28S:
    try:
        accession = seq.id.split('|')[0].split(':')[1]
    except IndexError:
        accession = seq.id.split('_')[0]
    accessions.append(accession)
accessions = sorted(set(accessions))
acc_to_name = {accession:reference.loc[reference['accession'] == accession,'name'].iloc[0] for accession in accessions}

#search for each name on open tree of life and retrieve information
ott_version = 'ott_version:' + str(requests.post('https://api.opentreeoflife.org/v3/taxonomy/about').json()['source'])
name_to_info =  {i:dict() for k,i in acc_to_name.iteritems()}

for k in name_to_info.iterkeys():
    query = ' '.join(k.split(' ')[0:2])
    r = requests.post('https://api.opentreeoflife.org/v3/tnrs/match_names',
                    data = {'names':[query, query], #since some names are replicated, sending unique names
                            'do_approximate_matching':True,
                            'context_name':'Arthropods'})
    #for some reason, OTT api strips 'sp.' from 'matched_name'. Below I try to adress it
    result = [j for j in r.json()['results'][0]['matches'] if (j['matched_name'] == query or
                                                                j['matched_name'] == query[:-4] or
                                                                j['matched_name'] == query[:-3])][0]
    name_to_info[k]['ott_id'] = str(result['taxon']['ott_id'])


for k,i in name_to_info.iteritems():
    r = requests.post('https://api.opentreeoflife.org/v3/taxonomy/taxon_info',
                data = {'ott_id': i['ott_id'], #id for taxon being searched
                        'include_lineage':True}) #include higher taxa
    for lineage in r.json()['lineage']:
        if lineage['rank'] == 'genus':
            name_to_info[k]['genus'] = lineage['name']
            name_to_info[k]['g_ott_id'] = str(lineage['ott_id'])
        if lineage['rank'] == 'subfamily':
            name_to_info[k]['subfamily'] = lineage['name']
        if lineage['rank'] == 'family':
            name_to_info[k]['family'] = lineage['name']
            name_to_info[k]['f_ott_id'] = str(lineage['ott_id'])
    if r.json()['rank'] == 'genus':
        name_to_info[k]['genus'] = r.json()['name']
        name_to_info[k]['g_ott_id'] = name_to_info[k]['ott_id']
    elif r.json()['rank'] == 'family':
        name_to_info[k]['family'] = r.json()['name']
        name_to_info[k]['f_ott_id'] = name_to_info[k]['ott_id']
    elif r.json()['rank'] == 'species':
        name_to_info[k]['s_ott_id'] = name_to_info[k]['ott_id']


#download misof tree and get taxonomic information for each tip
misof_tree = dendropy.Tree.get(url='https://api.opentreeoflife.org/v3/study/ot_211.nex?tip_label=ot:ottid', schema="nexus")

ott_ids_to_names = {i['ott_id']:k for k,i in name_to_info.iteritems()}
for leaf in misof_tree.leaf_node_iter():
    ott_id=leaf.taxon.label
    if ott_id in ott_ids_to_names.keys():
        taxname = ott_ids_to_names[ott_id]
        leaf.taxon.label = '|'.join(sorted([k + ':' + i for k,i in name_to_info[taxname].iteritems()]) + [ott_version,taxname])
        ott_ids_to_names['final_seq_name'] = leaf.taxon.label
    else:
        r = requests.post('https://api.opentreeoflife.org/v3/taxonomy/taxon_info',
            data = {'ott_id': ott_id, #id for taxon being searched
                    'include_lineage':True}) #include higher taxa
        tempdict = dict()
        tempdict['taxname'] = r.json()['name']
        for lineage in r.json()['lineage']:
            if lineage['rank'] == 'genus':
                tempdict['genus'] = lineage['name']
                tempdict['g_ott_id'] = str(lineage['ott_id'])
            if lineage['rank'] == 'family':
                tempdict['family'] = lineage['name']
                tempdict['f_ott_id'] = str(lineage['ott_id'])
        if r.json()['rank'] == 'genus':
            tempdict['genus'] = r.json()['name']
            tempdict['g_ott_id'] = ott_id
        leaf.taxon.label = '|'.join(sorted([k + ':' + i for k,i in tempdict.iteritems() if k != 'taxname']) + [ott_version,tempdict['taxname']])


#save full misof tree
misof_tree.write(path='misof_full_annotated.tre',schema='nexus')

#prune genera not found in sequences and save prunned version of tree
restart = True
while restart:
    for leaf in misof_tree.leaf_node_iter():
        try: #if no g_ott_id or g_ott_id not found in sequences, remove
            g_ott_id = [info.split(':')[1] for info in leaf.taxon.label.split('|')[:-1] if info.split(':')[0] == 'g_ott_id'][0]
        except IndexError:
            misof_tree.prune_taxa_with_labels([leaf.taxon.label])
            break
        if g_ott_id not in [i['g_ott_id'] for k,i in name_to_info.iteritems() if 'g_ott_id' in i.keys()]:
            misof_tree.prune_taxa_with_labels([leaf.taxon.label])
            break
    else:
        restart = False

#save prunned misof tree
misof_tree.purge_taxon_namespace()
misof_tree.write(path='misof_prun_annotated.tre',schema='nexus')


#loop through sequences, rename them and save
for seq in g18S:
    try:
        accession = seq.id.split('|')[0].split(':')[1]
    except IndexError:
        accession = seq.id.split('_')[0]
    seq.id = '|'.join(sorted([accession] + ['18S'] + [k + ':' + i for k,i in name_to_info[acc_to_name[accession]].iteritems()]) + [ott_version,acc_to_name[accession]])
    seq.name = ''
    seq.description = ''

for seq in g28S:
    try:
        accession = seq.id.split('|')[0].split(':')[1]
    except IndexError:
        accession = seq.id.split('_')[0]
    seq.id = '|'.join(sorted([accession] + ['28S'] + [k + ':' + i for k,i in name_to_info[acc_to_name[accession]].iteritems()]) + [ott_version,acc_to_name[accession]])
    seq.name = ''
    seq.description = ''
SeqIO.write(g18S, open('18S_renamed.fasta','w'),'fasta')
SeqIO.write(g28S, open('28S_renamed.fasta','w'),'fasta')


#now, produce versions of tree and alignment with genus name and genus ott id as the only identifiers, to be used by papara/raxml in taxonomic placement
for leaf in misof_tree.leaf_node_iter():
    tempdict = {field.split(':')[0]:field.split(':')[1] for field in leaf.taxon.label.split('|')[:-1]}
    leaf.taxon.label = tempdict['genus'] + '_' + tempdict['g_ott_id']

for seq in g18S:
    tempdict = {field.split(':')[0]:field.split(':')[1] for field in seq.id.split('|')[2:-1]}
    seq.id = tempdict['genus'] + '_' + tempdict['g_ott_id']

for seq in g28S:
    tempdict = {field.split(':')[0]:field.split(':')[1] for field in seq.id.split('|')[2:-1]}
    seq.id = tempdict['genus'] + '_' + tempdict['g_ott_id']

misof_tree.write(path='misof_papara.tre',schema='newick')
SeqIO.write(g18S, open('18S_papara.fasta','w'),'fasta')
SeqIO.write(g28S, open('28S_papara.fasta','w'),'fasta')
