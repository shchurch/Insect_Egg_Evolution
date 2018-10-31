#!/usr/bin/env python

### Created by Bruno de Medeiros (souzademedeiros@fas.harvard.edu), starting on 25-May-2017

### This script takes the Misof tree downloaded from open tree of life, already annotated with taxonomic information
### It returns a nexml tree with nodes named and estimated ages as annotations

import dendropy, pandas, sys, numpy
sys.path.insert(0,'../') #necessary to load parse_seqname
sys.path.insert(0,'../../') #necessary to load parse_seqname
from silva_to_alignment import parse_seqname


def annotate_tips(tip):
    annotations = parse_seqname(tip.taxon.label)
    try:
        tip.taxon.label = annotations['genus'] + '_' + annotations['g_ott_id']
    except KeyError:
        tip.taxon.label = annotations['family'] + '_' + annotations['f_ott_id']
    for i,x in annotations.iteritems():
        tip.annotations.add_new(i,x)
    
def build_tree_from_table(tree_as_table):
    tree = dendropy.Tree()
    tree.seed_node = dendropy.Node()
    
    #assuming the first line is the parent node
    tree.seed_node.new_child(label = tree_as_table.iloc[0,0])
    
    #now build the tree
    for i,row in tree_as_table.iterrows():
        current_node = tree.find_node_with_label(row.iloc[0])
        if row.iloc[1] not in set(tree_as_table['parent']):
            current_node.new_child(taxon = dendropy.Taxon(label = row.iloc[1]))
        else:
            current_node.new_child(label = row.iloc[1])
            
        
    tree.update_taxon_namespace()
    tree.update_bipartitions()
    
    return tree
    

if __name__ == "__main__":
    
    #first, read misof tree and rename tips
    misof_original = dendropy.Tree.get(path='misof_full_annotated.tre', schema = 'nexus', rooting='force-rooted')
    for tip in misof_original.leaf_node_iter():
        annotate_tips(tip)
        
    #now, find mrcas and add order names (not doing as taxon since some orders are represented by one tip only)
    misof_trans_table = pandas.DataFrame.from_csv('misof_name_translation.csv')
    for taxon_name in set(misof_trans_table['reference_taxon']):
        tips = [x['ott_genus'] + '_' + str(x['ott_id_genus']) for i,x in misof_trans_table.iterrows() if x['reference_taxon'] == taxon_name]
        mrca_node = misof_original.mrca(taxon_labels = set(tips))
        mrca_node.annotations.add_new('order', taxon_name)
        
    #load misof figure topology as a tree and replace names for the names in our tree
    tree_as_table = pandas.DataFrame.from_csv('misof_figure_tree_topology.tsv',sep='\t', index_col = None)
    misof_figure = build_tree_from_table(tree_as_table)
    for leaf in misof_figure.leaf_node_iter():
        if leaf.taxon.label in set(misof_trans_table['misof_figure_name']):
            leaf.taxon.label = [row['ott_genus'] + '_' + str(row['ott_id_genus']) for i,row in misof_trans_table.iterrows() if row['misof_figure_name'] == leaf.taxon.label][0]
        else:
            pass
    
    
    #now transfer node labels from misof figure to our tree
    for node in misof_original.postorder_internal_node_iter():
        tip_names = [leaf.taxon.label for leaf in node.leaf_iter()]
        try:
            figure_node = misof_figure.mrca(taxon_labels = tip_names)
            node.annotations.add_new('misof_node_number', figure_node.label)
        except KeyError:
            pass
        
    #finally, load time estimates from table and add annotations
    misof_time = pandas.DataFrame.from_csv('misof_File_s17.txt', sep = '\t')
    misof_time.index = misof_time.index.astype('string')
    
    
    misof_medians = dict()
    misof_2_5 = dict()
    misof_97_5 = dict()
    misof_all = dict()
    for i,row in misof_time.iterrows():
        misof_all[i] = row.dropna().as_matrix()
        misof_medians[i] = numpy.median(misof_all[i])
        misof_2_5[i] = numpy.percentile(misof_all[i], 2.5)
        misof_97_5[i] = numpy.percentile(misof_all[i], 97.5)
    
    for node in misof_original.postorder_internal_node_iter():
        ann_dict = node.annotations.values_as_dict()
        if 'misof_node_number' in ann_dict.keys() and ann_dict['misof_node_number'] in misof_time.index:
            node_number = ann_dict['misof_node_number']
            node.annotations.add_new('age_median', misof_medians[node_number])
            node.annotations.add_new('age_95%_interval', [misof_2_5[node_number],misof_97_5[node_number]])
            node.annotations.add_new('age_all', list(misof_all[node_number]))
            
    #save the fully annotated tree as nexml
    misof_original.update_taxon_namespace()
    misof_original.write(path='fully_annotated_misof.nexml',schema='nexml')
    
    #save table with taxa and estimated ages
    ages = []
    for node in misof_original.postorder_internal_node_iter():
        if node.taxon and node.taxon.label != 'Crustacea':
            ages.append({'taxa':node.taxon.label,
                        'median_age':'{:.2f}'.format(node.annotations.get_value('age_median')),
                        'age_95%_interval':'-'.join(['{:.2f}'.format(x) for x in node.annotations.get_value('age_95%_interval')])})
        elif 'misof_node_number' in node.annotations.values_as_dict().keys() and 'no_name' not in node.annotations.values_as_dict()['misof_node_number']:
            all_children = ','.join(sorted([x.taxon.label for x in node.postorder_internal_node_iter() if x.taxon is not None]))
            if all_children:
                all_previous_taxa = [age['taxa'] for age in ages]
                if all_children == all_previous_taxa: #since we are traversing in postorder, it will visit the mrca first
                    continue
                else:
                    ages.append({'taxa':all_children,
                                 'median_age':'{:.2f}'.format(node.annotations.get_value('age_median')),
                                 'age_95%_interval':'-'.join(['{:.2f}'.format(x) for x in node.annotations.get_value('age_95%_interval')])})
                    
    pandas.DataFrame(ages).to_csv('../ML_tree/calibrations.tsv', sep = '\t', index = False)
