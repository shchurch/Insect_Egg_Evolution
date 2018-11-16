#!/usr/bin/env python
#this is a simple script that converts fasta to phylip, keeping only accession as sequence identifier
#gene name as to be provided as command-line argument

from Bio import AlignIO
import sys
gene = sys.argv[1]

alignment = AlignIO.read(open(gene + '_epa_alignment.fasta','r'),'fasta')
for record in alignment:
	record.id = record.description.split('|')[0]
	record.description = ''
	record.name = ''
AlignIO.write(alignment,open(gene + '_epa_alignment.phy','w'),'phylip-relaxed')
