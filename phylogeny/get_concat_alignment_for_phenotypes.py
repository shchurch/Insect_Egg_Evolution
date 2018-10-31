#!/usr/bin/env python

### Created by Bruno de Medeiros (souzademedeiros@fas.harvard.edu), starting on 01-nov-2016
### The purpose of the script is to generate final phylogenetic datasets for genera and species, filtering out sequences with no phenotypic data
### The final matrix is written as a phylip file, for raxml
### It also makes tables with the overlap of genera and species between phenotypic and phylogenetic datasets
### This script is a little messy now, and needs optimization. In any case, it works.

### Note: this script is very messy, needs to be rewritten to accomodate more genes!


from Bio import SeqIO, AlignIO, Align
from pandas import DataFrame
from silva_to_alignment import parse_seqname
import dendropy, re, random

all_seqs_18S = [record for record in SeqIO.parse('18S/18S_aligned_final_trimmed.fasta','fasta')]
all_seqs_28S = [record for record in SeqIO.parse('28S/28S_aligned_final_trimmed.fasta','fasta')]

#mafft truncates seqnames, so we need to recover them from full sequences
#seqnames_18S = [record.description for record in SeqIO.parse('18S/18S_fullseqs_final.fasta','fasta')]
#seqnames_28S = [record.description for record in SeqIO.parse('28S/28S_fullseqs_final.fasta','fasta')]
#for record in all_seqs_18S:
#    record.description = [fullname for fullname in seqnames_18S if fullname.split('|')[0] == record.description.split('|')[0]][0]
#for record in all_seqs_28S:
#    record.description = [fullname for fullname in seqnames_28S if fullname.split('|')[0] == record.description.split('|')[0]][0]
    
#record genera and species in alignments    

genera_ott_in_alignments = set()
species_ott_in_alignments = set()

for record in all_seqs_18S:
    genera_ott_in_alignments.update({int(parse_seqname(record.description)['g_ott_id'])})
    try:
        species_ott_in_alignments.update({int(parse_seqname(record.description)['s_ott_id'])})
    except KeyError:
        pass
for record in all_seqs_28S:
    genera_ott_in_alignments.update({int(parse_seqname(record.description)['g_ott_id'])})
    try:
        species_ott_in_alignments.update({int(parse_seqname(record.description)['s_ott_id'])})
    except KeyError:
        pass
    
genera_ott_in_alignments = genera_ott_in_alignments - set([None])
species_ott_in_alignments = species_ott_in_alignments - set([None])

#make a species-level alignment
s_18S = []
for record in all_seqs_18S:
    if 's_ott_id' in record.description:
        s_18S.append(int(parse_seqname(record.description)['s_ott_id']))
    else:
        s_18S.append(None)

s_28S = []
for record in all_seqs_28S:
    if 's_ott_id' in record.description:
        s_28S.append(int(parse_seqname(record.description)['s_ott_id']))
    else:
        s_28S.append(None)



sp_seqs_18S = []
sp_seqs_28S = []

for sp in species_ott_in_alignments:
    seq_indexes_18S = [i for i in xrange(len(all_seqs_18S)) if s_18S[i] == sp]
    seq_indexes_28S = [i for i in xrange(len(all_seqs_28S)) if s_28S[i] == sp]
    if len(seq_indexes_18S) > 0:
        ngaps = [len(all_seqs_18S[i].seq) - len(re.findall('[ATCG]{1}',str(all_seqs_18S[i].seq))) for i in seq_indexes_18S]
        mingaps = min(ngaps)
        to_keep = random.choice([seq_indexes_18S[i] for i in xrange(len(seq_indexes_18S)) if ngaps[i] == mingaps])
        sp_seqs_18S.append(all_seqs_18S[to_keep][:])
    if len(seq_indexes_28S) > 0:
        ngaps = [len(all_seqs_28S[i].seq) - len(re.findall('[ATCG]{1}',str(all_seqs_28S[i].seq))) for i in seq_indexes_28S]
        mingaps = min(ngaps)
        to_keep = random.choice([seq_indexes_28S[i] for i in xrange(len(seq_indexes_28S)) if ngaps[i] == mingaps])
        sp_seqs_28S.append(all_seqs_28S[to_keep][:])




#make a genus-level alignment
g_18S = [int(parse_seqname(record.description)['g_ott_id']) for record in all_seqs_18S]
g_28S = [int(parse_seqname(record.description)['g_ott_id']) for record in all_seqs_28S]

seqs_18S = []
seqs_28S = []

for genus in genera_ott_in_alignments:
    seq_indexes_18S = [i for i in xrange(len(all_seqs_18S)) if g_18S[i] == genus]
    seq_indexes_28S = [i for i in xrange(len(all_seqs_28S)) if g_28S[i] == genus]
    if len(seq_indexes_18S) > 0:
        ngaps = [len(all_seqs_18S[i].seq) - len(re.findall('[ATCG]{1}',str(all_seqs_18S[i].seq))) for i in seq_indexes_18S]
        mingaps = min(ngaps)
        to_keep = random.choice([seq_indexes_18S[i] for i in xrange(len(seq_indexes_18S)) if ngaps[i] == mingaps])
        seqs_18S.append(all_seqs_18S[to_keep][:])
    if len(seq_indexes_28S) > 0:
        ngaps = [len(all_seqs_28S[i].seq) - len(re.findall('[ATCG]{1}',str(all_seqs_28S[i].seq))) for i in seq_indexes_28S]
        mingaps = min(ngaps)
        to_keep = random.choice([seq_indexes_28S[i] for i in xrange(len(seq_indexes_28S)) if ngaps[i] == mingaps])
        seqs_28S.append(all_seqs_28S[to_keep][:])


#make new alignments, only keeping genera and species for which we have phenotypic data
with open('../egg_database.txt', 'r') as infile:
    phenotypic_records = [eval(line) for line in infile]

genera_ottid = set()
species_ottid = set()


for record in phenotypic_records:
    try:
        genera_ottid.add(int(record['tax_cg_ott_id']))
    except (KeyError,ValueError):
        pass
    try:
        species_ottid.add(int(record['tax_cs_ott_id']))
    except (KeyError,ValueError):
        pass
genera_ottid = genera_ottid - set([None])
species_ottid = species_ottid - set([None])



genus2fam = dict()
for genus in genera_ottid:
    for record in phenotypic_records:
        if 'tax_cg_ott_id' in record.keys() and \
            record['tax_cg_ott_id'] != 'not_found' and \
            int(record['tax_cg_ott_id']) == genus and \
            'tax_family' in record.keys() and \
            record['tax_family'] is not None:
                
            genus2fam[genus] = record['tax_family']
            break
    else:
        genus2fam[genus] = None
        
sp2fam = dict()
for species in species_ottid:
    for record in phenotypic_records:
        if 'tax_cs_ott_id' in record.keys() and \
            record['tax_cs_ott_id'] != 'not_found' and \
            int(record['tax_cs_ott_id']) == species and \
            'tax_family' in record.keys() and \
            record['tax_family'] is not None:
                
            sp2fam[species] = record['tax_family']
            break
    else:
        sp2fam[species] = None
#print 'sp'        
#print species_ottid - set(sp2fam.keys())
#print 'genus'
#print genera_ottid - set(genus2fam.keys())


seqs_18S_with_data = []
for record in seqs_18S:
    if int(parse_seqname(record.description)['g_ott_id']) in genera_ottid:
        seqs_18S_with_data.append(record)

seqs_28S_with_data = []
for record in seqs_28S:
    if int(parse_seqname(record.description)['g_ott_id']) in genera_ottid:
        seqs_28S_with_data.append(record)
        
        
       
sp_seqs_18S_with_data = []
for record in sp_seqs_18S:
    if int(parse_seqname(record.description)['g_ott_id']) in genera_ottid:
        sp_seqs_18S_with_data.append(record)

sp_seqs_28S_with_data = []
for record in sp_seqs_28S:
    if int(parse_seqname(record.description)['g_ott_id']) in genera_ottid:
        sp_seqs_28S_with_data.append(record)



#write overlap tables
#for species
s_18S = []
for record in sp_seqs_18S:
    seqinfo = parse_seqname(record.description)
    if 's_ott_id' in seqinfo.keys():
        s_18S.append(int(parse_seqname(record.description)['s_ott_id']))
    else:
        s_18S.append(None)

s_28S = []
for record in sp_seqs_28S:
    seqinfo = parse_seqname(record.description)
    if 's_ott_id' in seqinfo.keys():
        s_28S.append(int(parse_seqname(record.description)['s_ott_id']))
    else:
        s_28S.append(None)

included_species_info = []
for sp in species_ottid:
    rec_18S = [sp_seqs_18S[i] for i in xrange(len(sp_seqs_18S)) if s_18S[i] == sp]
    rec_28S = [sp_seqs_28S[i] for i in xrange(len(sp_seqs_28S)) if s_28S[i] == sp]
    if len(rec_18S) > 1 or len(rec_28S) > 1:
        raise Exception('still more than one record per species')
    has_18S = len(rec_18S) > 0
    has_28S = len(rec_28S) > 0
    if has_18S or has_28S:
        try:
            seqinfo = parse_seqname(rec_18S[0].description)
        except:
            seqinfo = parse_seqname(rec_28S[0].description)
        try: #Zygentoma is considered subclass, with no order. Entognathes are at the class level
            order = seqinfo['order']
        except KeyError:
            try:
                order = seqinfo['subclass']
            except KeyError:
                order = seqinfo['class']
        try: #Sequences assembled de novo unfortunately lack the key current_name ->correct this later but for now use original_name
            name = seqinfo['current_name']
        except:
            name = seqinfo['original_name']
        included_species_info.append({'ott_order':order,'raxml_placement':seqinfo['raxml_placement'],'family':sp2fam[int(seqinfo['s_ott_id'])],'g_ott_id':genus,'genus':seqinfo['genus'],'full_name':name,'s_ott_id':seqinfo['s_ott_id'],'18S':has_18S,'28S':has_28S})
        

#for genera
included_genus_info = []
for genus in genera_ottid:
    rec_18S = [record for record in seqs_18S_with_data if int(parse_seqname(record.description)['g_ott_id']) == genus]
    rec_28S = [record for record in seqs_28S_with_data if int(parse_seqname(record.description)['g_ott_id']) == genus]
    if len(rec_18S) > 1 or len(rec_28S) > 1:
        raise Exception('still more than one record per genus')
    has_18S = len(rec_18S) > 0
    has_28S = len(rec_28S) > 0
    if has_18S or has_28S:
        try:
            seqinfo = parse_seqname(rec_18S[0].description)
        except:
            seqinfo = parse_seqname(rec_28S[0].description)
        try: #Zygentoma is considered subclass, with no order. Entognathes are at the class level
            order = seqinfo['order']
        except KeyError:
            try:
                order = seqinfo['subclass']
            except KeyError:
                order = seqinfo['class']

        included_genus_info.append({'ott_order':order,'raxml_placement':seqinfo['raxml_placement'],'family':genus2fam[int(genus)],'g_ott_id':genus,'genus':seqinfo['genus'],'18S':has_18S,'28S':has_28S})

#write tables to file
DataFrame.from_dict(included_species_info).to_csv(open('dataset_species_overlap.csv','w'),index = False, sep='\t')
DataFrame.from_dict(included_genus_info).to_csv(open('dataset_genus_overlap.csv','w'),index = False, sep='\t')

#write a concatenated genus alignment for raxml and partition file
all_labels = set()

for record in seqs_18S_with_data:
    seqinfo = parse_seqname(record.description)
    record.id = seqinfo['genus'] + '_' + seqinfo['g_ott_id']
    all_labels.update({record.id})
    record.description = ''
    record.name = ''

for record in seqs_28S_with_data:
    seqinfo = parse_seqname(record.description)
    record.id = seqinfo['genus'] + '_' + seqinfo['g_ott_id']
    all_labels.update({record.id})
    record.description = ''
    record.name = ''
    

#now remove column with gaps only
SeqIO.write(seqs_18S_with_data,open('18S/genus_alignment_18S.fasta','w'),'fasta')
SeqIO.write(seqs_28S_with_data,open('28S/genus_alignment_28S.fasta','w'),'fasta')

seqs_18S_with_data = AlignIO.read(open('18S/genus_alignment_18S.fasta','r'),'fasta')
seqs_28S_with_data = AlignIO.read(open('28S/genus_alignment_28S.fasta','r'),'fasta')


i = 0
temp_align = seqs_18S_with_data[:,0:0]
for j in xrange(seqs_18S_with_data.get_alignment_length()):
    if set(seqs_18S_with_data[:,j]) == set('-'):
        if j>i+1:
            temp_align = temp_align + seqs_18S_with_data[:,i+1:j]
        i = j
seqs_18S_with_data = temp_align

i = 0
temp_align = seqs_28S_with_data[:,0:0]
for j in xrange(seqs_28S_with_data.get_alignment_length()):
    if set(seqs_28S_with_data[:,j]) == set('-'):
        if j>i+1:
            temp_align = temp_align + seqs_28S_with_data[:,i+1:j]
        i = j
seqs_28S_with_data = temp_align
    
AlignIO.write(seqs_18S_with_data,open('18S/genus_alignment_18S.fasta','w'),'fasta')
AlignIO.write(seqs_28S_with_data,open('28S/genus_alignment_28S.fasta','w'),'fasta')

#now use dendropy to concatenate
all_labels = sorted(all_labels)
all_taxa = dendropy.TaxonNamespace(all_labels,label='taxa')

matrix_18S = dendropy.DnaCharacterMatrix.get(path='18S/genus_alignment_18S.fasta',schema='fasta',taxon_namespace=all_taxa)
matrix_18S.pack(matrix_18S.state_alphabets[0].gap)
matrix_28S = dendropy.DnaCharacterMatrix.get(path='28S/genus_alignment_28S.fasta',schema='fasta',taxon_namespace=all_taxa)
matrix_28S.pack(matrix_28S.state_alphabets[0].gap)


#write a concatenated genus alignment for raxml and partition file
dendropy.DnaCharacterMatrix.concatenate([matrix_18S,matrix_28S]).write(path='ML_tree/genus_concatenated.phy', schema='phylip') #using fasta for now because phylip results in some bug when reading by dendropy

with open('./ML_tree/genus_partitions','w') as outfile:
    outfile.write('DNA, 18S=1-' + str(matrix_18S.sequence_size) + '\n')
    outfile.write('DNA, 28S=' + str(matrix_18S.sequence_size + 1) + '-' + str(matrix_18S.sequence_size + matrix_28S.sequence_size))







#finally, repeat everything, building a species dataset

#write a concatenated species alignment for raxml and partition file
all_labels = set()

def get_species_name(original_name):
    return '_'.join(original_name.split('(')[0].strip().split(' ')[0:2])

for record in sp_seqs_18S_with_data:
    seqinfo = parse_seqname(record.description)
    record.id = get_species_name(seqinfo['original_name']) + '_' + seqinfo['s_ott_id']    
    all_labels.update({record.id})
    record.description = ''
    record.name = ''

for record in sp_seqs_28S_with_data:
    seqinfo = parse_seqname(record.description)
    record.id = get_species_name(seqinfo['original_name']) + '_' + seqinfo['s_ott_id']   
    all_labels.update({record.id})
    record.description = ''
    record.name = ''

#now remove columns with gaps only
SeqIO.write(sp_seqs_18S_with_data,open('18S/species_alignment_18S.fasta','w'),'fasta')
SeqIO.write(sp_seqs_28S_with_data,open('28S/species_alignment_28S.fasta','w'),'fasta')

sp_seqs_18S_with_data = AlignIO.read(open('18S/species_alignment_18S.fasta','r'),'fasta')
sp_seqs_28S_with_data = AlignIO.read(open('28S/species_alignment_28S.fasta','r'),'fasta')

i = 0
temp_align = sp_seqs_18S_with_data[:,0:0]
for j in xrange(sp_seqs_18S_with_data.get_alignment_length()):
    if set(sp_seqs_18S_with_data[:,j]) == set('-'):
        if j>i+1:
            temp_align = temp_align + sp_seqs_18S_with_data[:,i+1:j]
        i = j
sp_seqs_18S_with_data = temp_align
        

i = 0
temp_align = sp_seqs_28S_with_data[:,0:0]
for j in xrange(sp_seqs_28S_with_data.get_alignment_length()):
    if set(sp_seqs_28S_with_data[:,j]) == set('-'):
        if j>i+1:
            temp_align = temp_align + sp_seqs_28S_with_data[:,i+1:j]
        i = j
sp_seqs_28S_with_data = temp_align
    
AlignIO.write(sp_seqs_18S_with_data,open('18S/species_alignment_18S.fasta','w'),'fasta')
AlignIO.write(sp_seqs_28S_with_data,open('28S/species_alignment_28S.fasta','w'),'fasta')

#now use dendropy to concatenate
all_labels = sorted(all_labels)
all_taxa = dendropy.TaxonNamespace(all_labels,label='taxa')

matrix_18S = dendropy.DnaCharacterMatrix.get(path='18S/species_alignment_18S.fasta',schema='fasta',taxon_namespace=all_taxa)
matrix_18S.pack(matrix_18S.state_alphabets[0].gap)
matrix_28S = dendropy.DnaCharacterMatrix.get(path='28S/species_alignment_28S.fasta',schema='fasta',taxon_namespace=all_taxa)
matrix_28S.pack(matrix_28S.state_alphabets[0].gap)

#now use biopython to remove alignment columns that are gaps only
dendropy.DnaCharacterMatrix.concatenate([matrix_18S,matrix_28S]).write(path='ML_tree/species_concatenated.phy', schema='phylip') #using fasta for now because phylip results in some bug when reading by dendropy

with open('./ML_tree/species_partitions','w') as outfile:
    outfile.write('DNA, 18S=1-' + str(matrix_18S.sequence_size) + '\n')
    outfile.write('DNA, 28S=' + str(matrix_18S.sequence_size + 1) + '-' + str(matrix_18S.sequence_size + matrix_28S.sequence_size))

