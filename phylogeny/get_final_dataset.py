#!/usr/bin/env python

### Created by Bruno de Medeiros (souzademedeiros@fas.harvard.edu), starting on 01-nov-2016
### The purpose of the script is to generate a final phylogenetic dataset, filtering out sequences with mismatches between taxonomy and phylogenetic placement
### Finally, genera included are tabulated

from Bio import SeqIO
from pandas import DataFrame
from collections import defaultdict
from silva_to_alignment import parse_seqname
import argparse, dendropy, json, numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('gene', help = '18S or 28S')
parser.add_argument('-i','--info', default = './raxml_epa/ref_taxonomy.csv', help = 'path to table with taxonomic information about reference sequences')

#args = parser.parse_args(['28S'])
args = parser.parse_args()




#1 - get an annotated reference tree
raxml = json.load(open('./raxml_epa/RAxML_portableTree.'+ args.gene + '.jplace', 'r'))
tree = dendropy.Tree.get(data = '[&R] ' + raxml['tree'], schema = 'newick', preserve_underscores=True, is_parse_jplace_tokens = True)
info = DataFrame.from_csv(args.info, index_col=None)

clades = dict()
for clade in set(info['reference_taxon']):
    tips = set(info['reference_name'].loc[info['reference_taxon'] == clade]) - set([np.nan])
    clades.update({clade:tips})

for clade, tips in clades.iteritems():
    mrca_node = tree.mrca(taxon_labels=tips)
    mrca_node.annotations.add_new('clade',clade)

tree.seed_node.annotations.add_new('clade','Hexapoda')

for node in tree.postorder_node_iter():
    if 'clade' not in node.annotations.values_as_dict().keys():
        node.annotations.add_new('clade',None)

#2 - load query sequences
## loop through query sequences and add information about taxonomic placement
## if taxonomic placement matches ott order, put on list matches. If not, on mismatches

misof_to_ott = defaultdict(set)
for index in info.index:
    misof_to_ott[info.loc[index,'reference_taxon']].update([info.loc[index,'ott_order']])


matches = []
mismatches = []

original_seqs = dict()
#reading the unaligned sequences because we will do another alignment
for record in SeqIO.parse(args.gene + '/query_unaligned.fasta','fasta'):
    seq_info = parse_seqname(record.description)
    original_seqs.update({seq_info['accession']:{'record':record, 'info':seq_info}})

placed_seqs = dict()
for record in SeqIO.parse('raxml_epa/' + args.gene + '_epa_alignment.fasta','fasta'):
    if '|' in record.description: #only queries will have this symbol
        seq_info = original_seqs[record.description.split('|')[0]]['info']
        full_seq = original_seqs[record.description.split('|')[0]]['record']
        placed_seqs.update({seq_info['accession']:{'record':full_seq,'info':seq_info}})


for placement in raxml['placements']:
    if len(placement['n']) > 1:
        raise Exception('more than onse sequence in this placement, check')
    p_cum = defaultdict(float)
    for p_info in placement['p']:
        edge = tree.edge_index[p_info[0]]
        node = edge.head_node
        while True:
            label = node.annotations.get_value('clade')
            if label in clades.keys() + ['Hexapoda']:
                p_cum[label] += p_info[2]
                break
            else:
               node = node.parent_node
               continue

    probs = [p for k,p in p_cum.iteritems()]
    maxprobs = max(probs)
    best_placement = [k for k,p in p_cum.iteritems() if p == maxprobs][0]
    seqinfo = placed_seqs[placement['n'][0]]['info']
    seqinfo.update({'raxml_placement':best_placement})
    seq = placed_seqs[placement['n'][0]]['record']
    seq.description = ''
    seq.name = ''
    original_name = seqinfo.pop('original_name')
    accession = seqinfo.pop('accession')
    ott_version = 'ott_version:' + seqinfo.pop('ott_version')
    seq.id = '|'.join([accession] + sorted([k + ':' + str(x) for k,x in seqinfo.iteritems()]) + [ott_version,original_name])

    #since our concept of order is not consistently classified as order in ott, look at order, subclass and class fields for a match in raxml_placement
    if best_placement != 'Hexapoda' and 'order' in seqinfo.keys() and seqinfo['order'] in misof_to_ott[best_placement]:
        matches.append(seq)
    elif best_placement != 'Hexapoda' and 'subclass' in seqinfo.keys() and seqinfo['subclass'] in misof_to_ott[best_placement]:
        matches.append(seq)
    elif best_placement != 'Hexapoda' and 'class' in seqinfo.keys() and seqinfo['class'] in misof_to_ott[best_placement]:
        matches.append(seq)
    else:
        mismatches.append(seq)



SeqIO.write(matches,open(args.gene + '/' + args.gene + '_fullseqs_final.fasta','w'),'fasta')
SeqIO.write(mismatches,open(args.gene + '/' + args.gene + '_rejected.fasta','w'),'fasta')
