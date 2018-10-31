#!/usr/bin/env python

### Created by Bruno de Medeiros (souzademedeiros@fas.harvard.edu), starting on Nov-2016
### Because some 18S and 28S sequences were de novo assembled with no strand information, some of them had to be reverse-complemented for alignment.
### Since this was done in Geneious with alignments for papara, with genus name only as sequence names, I wrote this script to go back to the original sequence file
### and replace sequences that were reverse-complemented


from Bio import SeqIO
import argparse


#first, read input and save information in a pandas dataframe
parser = argparse.ArgumentParser()
parser.add_argument('-1', '--old18S', help = 'path to input 18S fasta file with full names', default = '18S_renamed.fasta')
parser.add_argument('-2', '--old28S', help = 'path to input 28S fasta file with full names', default = '28S_renamed.fasta')
parser.add_argument('-a', '--new18S', help = 'path to input 18S fasta file in correct orientation', default = '18S_papara.fasta')
parser.add_argument('-b', '--new28S', help = 'path to input 28S fasta file in correct orientation', default = '28S_papara.fasta')

args = parser.parse_args()
#args = parser.parse_args([])

#read sequence files and reference table
old18S_names = [record.description for record in SeqIO.parse(args.old18S,'fasta')]
old28S_names = [record.description for record in SeqIO.parse(args.old28S,'fasta')]
new18S = [record for record in SeqIO.parse(args.new18S,'fasta')]
new28S = [record for record in SeqIO.parse(args.new28S,'fasta')]

#open old sequences and create a dictionary of sequence names keyed by genus ott id
dict_18S = {{field.split(':')[0]:field.split(':')[1] for field in name.split('|')[2:-1]}['g_ott_id']:name for name in old18S_names}
dict_28S = {{field.split(':')[0]:field.split(':')[1] for field in name.split('|')[2:-1]}['g_ott_id']:name for name in old28S_names}

#loop through papara sequences and replace short name with full name
for seq in new18S:
    seq.id = dict_18S[seq.id.split('_')[1]]
    seq.name = ''
    seq.description = ''

for seq in new28S:
    seq.id = dict_28S[seq.id.split('_')[1]]
    seq.name = ''
    seq.description = ''

SeqIO.write(new18S, open('18S_renamed.fasta','w'),'fasta')
SeqIO.write(new28S, open('28S_renamed.fasta','w'),'fasta')
