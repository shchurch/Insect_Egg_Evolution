#!/usr/bin/env python

### Created by Bruno de Medeiros (souzademedeiros@fas.harvard.edu), starting on 03-jun-17

### This script trims an alignment to the region with a given minimum coverage of ingroup taxa
### It then removes sequences with more than 95% missing or ambiguous

import argparse, pandas, sys, dendropy
from collections import Counter
from os.path import basename, dirname

#after joining some sequences from misplaced, we have again multiple sequences per taxonomic id
#this function finds them and keeps only the longest
def find_region_of_ingroup_high_coverage(alignment, outgroup_labels, minimum_coverage = 0.3, minimum_w_data=10):
    alignment = alignment.clone(2)
    out_taxa = [taxon for taxon in alignment.taxon_namespace if any([outgroup in taxon.label for outgroup in outgroup_labels])]
    alignment.remove_sequences(out_taxa)
    gapamb = ['-','N']

    #find alignment start and end
    idx_above_mincov = []
    idx_above_min_data = []
    nseqs = len(alignment.sequences())
    for i in xrange(alignment.sequence_size):
        temp_gapamb = 0
        for seq in alignment.sequences():
            if str(seq[i]) in gapamb:
                temp_gapamb += 1

        if 1 - float(temp_gapamb)/nseqs >= minimum_coverage:
            idx_above_mincov.append(i)
        if nseqs - temp_gapamb >= minimum_w_data:
            idx_above_min_data.append(i)

    start = min(idx_above_mincov)
    end = max(idx_above_mincov)
    idx = [i for i in idx_above_min_data if i >= start and i <= end]


    print str(idx)
    return idx
    
def number_not_missing(seq):
    counts = Counter(seq)
    not_missing = sum([counts[nuc] for nuc in ['A','C','T','G','R','Y','M','S','W','K']])
    return not_missing
    

def find_taxa_with_little_data(alignment, min_data = 100):
    taxa_to_drop = list()
    for taxon in alignment.taxon_namespace:
        seq = alignment[taxon]
        if number_not_missing(seq.symbols_as_string()) < min_data:
            taxa_to_drop.append(taxon)
    return taxa_to_drop
            

    
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-a','--alignment', help = 'alignment in fasta format')
    parser.add_argument('-o','--outgroups', help = 'table with outgroup taxonomic ids for each clade')
    parser.add_argument('-c','--coverage', help = 'minimum coverage on the edge of alignments', default = 0.2)
    parser.add_argument('-m','--min-data', help = 'minimum number of non-missing nucleotides to keep a sample', default = 100)
    parser.add_argument('-s','--min-seqs', help  = 'minimum number of non-gap sequences to keep a position in the middle of the alignment', default = 10)


    args = parser.parse_args()
    #args = parser.parse_args([])
    
    print 'Trimming ' + args.alignment
    
    if args.outgroups:
        with open(args.outgroups,'r') as infile:
            as_table = pandas.read_csv(infile)
            all_outgroups = {x['higher_taxon']:str(x['ncbi_id']) for i,x in as_table.iterrows()}
    else:
        all_outgroups = dict()

    outgroup_labels = [x for i,x in all_outgroups.iteritems() if i not in args.alignment]
    
    alignment = dendropy.DnaCharacterMatrix.get(path = args.alignment, schema = 'fasta')
    
    
    #trim alingnment to region with high coverage
    #remove all positions with less than 10 terminals with data
    idx = find_region_of_ingroup_high_coverage(alignment, 
                                               outgroup_labels,
                                               minimum_coverage = float(args.coverage),
                                               minimum_w_data = int(args.min_seqs))
    
    new_alignment = alignment.export_character_indices(idx)

    #now remove species with less than 100 bp non-missing remaining
    taxa_to_drop = find_taxa_with_little_data(new_alignment, min_data = float(args.min_data))

    if taxa_to_drop:
        print args.alignment
        print 'The following taxa were removed due to excess missing data:'
        print '\n'.join([taxon.label for taxon in taxa_to_drop])
        new_alignment.remove_sequences(taxa_to_drop)
    
    #write output
        
    outpath = './' + dirname(args.alignment) + '/'  + basename(args.alignment).split('.')[0] + '_trimmed.fasta' 
    with open(outpath, 'w') as outfile:
        new_alignment.write(file = outfile, schema = 'fasta')

