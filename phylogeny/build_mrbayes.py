#!/usr/bin/env python

#This script compares Rainford et al tree to our dataset and generates a family-level consensus tree


import dendropy, pandas, requests, sys, argparse, random, os
import numpy as np
sys.path.insert(0,'../') #necessary to load otl_tnrs
from get_taxonomy import otl_tnrs, otl_taxon

#function updates taxon names in rainford according to open tree of life, and returns None for non-monophyletic groups
#monophyly is assumed if a node exists in open tree of life
def update_taxon(taxon, context = 'Arthopods'):
    r = otl_tnrs(taxon, context)
    if r.json()['results']: #if results found, return the best
        results = r.json()['results'][0]['matches']
        scores = [results[i]['score'] for i in range(len(results))] #make a list with matches' scores
        best = scores.index(max(scores)) #returns index for result with highest score. If more than one, keeps first
        if max(scores) == 1 and not results[best]['is_synonym']:
            best_name = results[best]['taxon']['name']
            best_id = results[best]['taxon']['ott_id']
            best_res = results[best]
        else: #if no perfect match, try to look synonyms
            syns = [res for res in results if taxon in res['taxon']['synonyms']]
            scores = [results[i]['score'] for i in range(len(syns))]
            best = scores.index(max(scores))
            best_res = syns[best]
            best_name = syns[best]['taxon']['name']
            best_id = syns[best]['taxon']['ott_id']
        #check if match was to family
        #if not, check if it is now a taxon within a family
        #if not, return None
        if best_res['taxon']['rank'] != 'family':
            r = otl_taxon(best_id)
            for lineage in r.json()['lineage']:
                if lineage['rank'] == 'family':
                    best_name = lineage['name']
                    best_id = str(lineage['ott_id'])
            else:
                return None

        r1 = requests.post('https://api.opentreeoflife.org/v3/tree_of_life/node_info',
                        data = {'ott_id':best_id, #id for taxon being searched
                                'include_lineage':False}) #include higher taxa
        if 'node_id' in r1.json().keys():
            return best_name
        else:
            return None

#this function removes columns with only gaps for a dendropy charmatrix  
#it preserves the original character subsets                      
def remove_gaps(charmatrix):
    idx_to_keep = []
    gapamb = ['-','N']
    for i in xrange(charmatrix.sequence_size):
        for seq in charmatrix.sequences():
            if str(seq[i]) not in gapamb:
                idx_to_keep.append(i)
                break
    newmatrix = charmatrix.export_character_indices(idx_to_keep)
    #now transfer subsets
    new_subsets = {name:{'last':0,'N':0} for name in charmatrix.character_subsets.viewkeys()}
    for idx in idx_to_keep:
        for name in charmatrix.character_subsets.viewkeys():
            if idx in charmatrix.character_subsets[name].character_indices:
                new_subsets[name]['last'] = max([new_subsets[name]['last'], idx])
                new_subsets[name]['N'] += 1
                break
    ordered_subsets = {x['last']:{'name':i, 'N':x['N']} for i, x in new_subsets.iteritems()}
    first=0
    for last in sorted(ordered_subsets.keys()):
        new_last = first + ordered_subsets[last]['N']
        newmatrix.new_character_subset(label = ordered_subsets[last]['name'], 
                                       character_indices = range(first,new_last))
        first = new_last
    
    return newmatrix

def make_big_tree(alignment, align_data, order_tree, subordinal_trees, outfolder, palaeoptera_paraphyletic = False):
    print 'Making mrbayes file for big tree'
    
    if not os.path.exists(outfolder + '/'):
        os.makedirs(outfolder)
    
    with open(outfolder + '/mrbayes_input.nex','w') as mrbayes_out:
        #first, write nexus without mrbayes block
        nexus_string = alignment.as_string(schema = 'nexus', unquoted_underscores=True)
        nexus_splits = nexus_string.split('BEGIN SETS;')
        mrbayes_out.write(nexus_splits[0])
        
        #now, start mrbayes block and add character partitions with GTR + GAMMA + I model
        mrbayes_out.write('\n\n' + 'BEGIN MRBAYES;' +  '\n\n')
        mrbayes_out.write(nexus_splits[1].split('END')[0] + '\n')
        mrbayes_out.write('partition by_gene = ' + 
                        str(len(alignment.character_subsets)) + 
                        ':' + 
                        ', '.join(alignment.character_subsets.keys()) +
                        ';\n')
        mrbayes_out.write('set partition = by_gene;\n')
        mrbayes_out.write('lset applyto=(all) nst=6 rates=invgamma;\n')
        mrbayes_out.write('unlink statefreq=(all) revmat=(all) shape=(all) pinvar=(all);\n')
        mrbayes_out.write('prset ratepr = variable;\n\n')
        
        #now add taxon constraints and calibrations
        #these lists will record taxset, constraint and calibrate stataments, to order them in the output
        tax_statements = []
        con_statements = []
        cal_statements = []
                    
        #hard constraints on all monophyletic families
#        families = set([str(x) for x in align_data.family if x is not np.nan])
#        align_name = {genus:genus + '_' + str(align_data.loc[[x == genus for x in align_data.genus],'g_ott_id'].iloc[0]) for genus in set(align_data.genus)}
        
        constraint_idx = 1
        for order,tree in subordinal_trees.iteritems():
           
            order_children_tips = set([x['genus'] + '_' + str(x['g_ott_id']) for i,x in align_data.iterrows() if x['raxml_placement'] == order])
            
                        
            try:
                order_children_families = [leaf.taxon.label for leaf in tree.leaf_nodes()]
            except:
                order_children_families = []
            
                
            #tax statements and constraints for families
            #constrain the current order
            tax_statements.append('taxset ' + order + ' = ' + ' '.join(order_children_tips))
            con_statements.append('constraint ' + order + ' hard = ' + order)
            
                            
            #family-level constraints only if families represented in rainford tree        
            if order_children_families:  
                #constrain monophyletic families
                for family in set([x['family'] for i,x in align_data.iterrows() if x['raxml_placement'] == order and x['family'] is not np.nan]):
                    family_tips = set([x['genus'] + '_' + str(x['g_ott_id']) for i,x in align_data.iterrows() if x['family'] == family])
                    tax_statements.append('taxset ' + family + ' = ' + ' '.join(family_tips))
                    if len(family_tips) > 1:
                        if family in rainford_reverse.keys() or update_taxon(family) is not None:#if family in Rainford dictionary, it was already tested as monophyletic
                            con_statements.append('constraint ' + family + ' hard = ' +  ' '.join(family_tips))
            
                #partial constraints on inter-familial relationships
                for node in tree.postorder_internal_node_iter():
                    #positive part of constraint
                    children_families = set([leaf.taxon.label for leaf in node.leaf_nodes()])
                    
                    #negative part of constraint
                    families_not = [leaf.taxon.label for leaf in tree.leaf_nodes() if leaf.taxon.label not in children_families]
                    
                    if families_not:
                        con_statement = 'constraint ' + 'cons_' + str(constraint_idx) + ' partial = ' + \
                                        ' '.join(sorted(children_families)) + \
                                        ' : ' + \
                                        ' '.join(sorted(families_not))
                        constraint_idx += 1
                        con_statements.append(con_statement)
        
        # for family in families:    
        #     genera = set(align_data.loc[[x == family for x in align_data.family],'genus'])
        #     genera = set([x['genus'] + '_' + str(x['g_ott_id']) for i,x in align_data.iterrows() if x['family'] == family]) 
        #     
        #     set([align_name[genus] for genus in genera])
        #     tax_statement = 'taxset ' + family + ' = ' + ' '.join(genera)
        #     con_statement = 'constraint ' + family + ' hard = ' + family
        #     #write taxset whether or not monophyletic
        #     tax_statements.append(tax_statement)
        # 
        #     #constraint only makes sense if more than one tip    
        #     if len(genera) > 1:        
        #     #check if monophyletic
        #         if update_taxon(family) is not None:#if None, family is not monophyletic in open tree of life
        #             con_statements.append(con_statement)
        
        
        # #create a taxset for each order, write hard constraint
        # #create taxsets and constraints for interfamilial relationships
        # constraint_idx = 1
        # for order,tree in subordinal_trees.iteritems():
        #     #hard constraint on orders, create a taxset for each order
        #     order_children_tips = set([x['genus'] + '_' + str(x['g_ott_id']) for i,x in align_data.iterrows() if x['raxml_placement'] == order])
        #     tax_statements.append('taxset ' + order + ' = ' + ' '.join(order_children_tips))
        #     con_statement = 'constraint ' + order + ' hard = ' + order
        #     
        #     #partial constraints on inter-familial relationships
        #     if tree:
        #         for node in tree.postorder_internal_node_iter():
        #             #positive part of constraint
        #             children_families = set([leaf.taxon.label for leaf in node.leaf_nodes()])
        #             
        #             #negative part of constraint
        #             families_not = [leaf.taxon.label for leaf in tree.leaf_nodes() if leaf.taxon.label not in children_families]
        #             if families_not:
        #                 con_statement = 'constraint ' + 'cons_' + str(constraint_idx) + ' partial = ' + \
        #                                 ' '.join(sorted(children_families)) + \
        #                                 ' : ' + \
        #                                 ' '.join(sorted(families_not))
        #                 constraint_idx += 1
        #                 con_statements.append(con_statement)
    
        
                    
        #hard constraints inter-ordinal relationships
        #calibrations for inter-ordinal relationships
        constraint_idx = 1
        for node in order_tree.preorder_internal_node_iter():
            children_orders = [leaf.taxon.label for leaf in node.leaf_nodes()]
            
            if set(children_orders) == set(['Ephemeroptera', 'Odonata']) and palaeoptera_paraphyletic is True:
                #if palaeoptera is to be treated as paraphyletic, only add constraints and ignore calibrations
                orders_external = ['Collembola', 'Protura', 'Diplura', 'Zygentoma', 'Archaeognatha']                    
                
                #constraint for Odonata + Neoptera
                con_statement = 'constraint ' + 'order_cons_' + str(constraint_idx) + ' hard = ' + ' '.join(sorted(set(subordinal_trees.keys()) - set(orders_external) - set(['Ephemeroptera'])))
                constraint_idx += 1
                #constraint for Ephemeroptera + Odonata + Neoptera
                con_statement = 'constraint ' + 'order_cons_' + str(constraint_idx) + ' hard = ' + ' '.join(sorted(set(subordinal_trees.keys()) - set(orders_external)))
                constraint_idx += 1
                continue
                
                
            else:
                con_statement = 'constraint ' + 'order_cons_' + str(constraint_idx) + ' hard = ' + ' '.join(sorted(children_orders))
            
            con_statements.append(con_statement)
            
            range_ages = []
            if order_tree.annotations.get_value('source') == 'rainford':
                #to get calibration statements, calibrations have to be searched in the full rainford tree
                children_tips = set()
                for order in children_orders:
                    if subordinal_trees[order] is None:
                        children_tips.add(order)
                    else:
                        children_tips = children_tips | set([leaf.taxon.label for leaf in subordinal_trees[order].leaf_nodes()])
                        
                rainford_names = [rainford_reverse[family] for family in children_tips]
                node_rainford = rainford_tree.mrca(taxon_labels=rainford_names)
                range_ages = [float(x) for x in node_rainford.annotations.values_as_dict()['age_95%HPD']]
            
                
            elif order_tree.annotations.get_value('source') == 'misof':
                #for the misof tree, calibrations are already there as annotations
                try:
                    range_ages = [float(x) for x in node.annotations.get_value('age_95%_interval').split()]
                except AttributeError:
                    range_ages = []
                
            if range_ages:   
                cal_min = min(range_ages)
                cal_max = max(range_ages)  
                cal_statement = 'calibrate ' + 'order_cons_' + str(constraint_idx) + ' = ' + \
                                'uniform(' + '{:.2f}'.format(cal_min) +  ',' + '{:.2f}'.format(cal_max) + ')'
                cal_statements.append(cal_statement)
            else:
                print 'No calibration for the ancestor of: ' + ','.join(children_orders)
            constraint_idx += 1 
        
            
        #finally, write all taxset, constraint and calibrate statements  
        #add tree and calibration priors

        mrbayes_out.write('prset brlenspr = clock:birthdeath;\n')
        constraints = [x.split(' ')[1] for x in con_statements]
        mrbayes_out.write(';\n'.join(tax_statements) + ';\n')
        mrbayes_out.write('\n\n')
        mrbayes_out.write(';\n'.join(con_statements)+ ';\n')
        mrbayes_out.write('\n\n')
        mrbayes_out.write(';\n'.join(cal_statements)+ ';\n')
        mrbayes_out.write('\n\n')
        mrbayes_out.write('prset topologypr = constraints(' + ','.join(constraints) + ');\n') 
        mrbayes_out.write('link topology=(all) brlens=(all);\n\n')
        
        
        priors = list()
        
        priors.append('sampleprob = 0.001') #we sampled about 1000 insect species form about 1M
        priors.append('speciationpr = exp(13.5)') #from the model with no shifts in condamine et al (2016) scientific reports (https://www.nature.com/articles/srep19208)
        priors.append('extinction = beta(2,2000)') #this puts extinction rate on the same order of magnitude as inferred in their models
        priors.append('nodeagepr = calibrated')
        priors.append('clockratepr = exponential(1000)') 
        priors.append('clockvarpr = igr') 
        priors.append('igrvarpr = exp(50)') 
        
                
        mrbayes_out.write('prset ' + ' '.join(priors) + ';\n')
        mrbayes_out.write('\n\n')
                
        #finally, write mcmc command
        mrbayes_out.write('showmodel;\n')
        mrbayes_out.write('mcmcp ngen=100000000 diagnfreq=10000 relburnin=yes burninfrac=0.25 printfreq=50000 samplefreq=10000 checkfreq=100000 savebrlens=yes;\n')
        mrbayes_out.write('[mcmcp append=yes;]\n') #uncomment later in the mrbayes code if restart is needed
        mrbayes_out.write('mcmc; sumt; sump; \n\n END;')
        

def make_orders_trees(alignment, align_data, order_tree, subordinal_trees, outfolder, palaeoptera_paraphyletic = False):

    if not os.path.exists(outfolder + '/'):
        os.makedirs(outfolder)

    #and several mrbayes files, one per order:
    print 'Making mrbayes files for each order'  
    #first, select randomly one representative per order, to be used as outgroups
    print 'Selecting outgroups:'  
    outgroups = dict()
    for order, tree in subordinal_trees.iteritems():
        order_children_tips = [x['genus'] + '_' + str(x['g_ott_id']) for i,x in align_data.iterrows() if x['raxml_placement'] == order and x['18S'] and x['28S']]
        if not order_children_tips:
            print order + ': no tip with both 18S and 28S to be used as outgroup, choosing one with 18S only.'
            order_children_tips = [x['genus'] + '_' + str(x['g_ott_id']) for i,x in align_data.iterrows() if x['raxml_placement'] == order and x['18S']]
        
        random.seed(order) #this will make sure the same outgroups are always chosen
        outgroups[order] = random.choice(order_children_tips)
        print order + ': ' + outgroups[order]
    
    
    #now, make a nexus file per order
    for order,tree in subordinal_trees.iteritems():
        
        order_children_tips = set([x['genus'] + '_' + str(x['g_ott_id']) for i,x in align_data.iterrows() if x['raxml_placement'] == order])
        print 'Making mrbayes file for ' + order 
        tax_statements = []
        con_statements = []
        cal_statements = []

        
        if len(order_children_tips) < 2:
            print order + ' has less than two tips in our dataset, skipping'
        
        
        else:
            with open(outfolder + '/' + order + '.nex', 'w') as mrbayes_out:
                #hard constraint on orders, create a taxset for each order            
                try:
                    order_children_families = [leaf.taxon.label for leaf in tree.leaf_nodes()]
                except:
                    order_children_families = []
                    
                
                
                tips_in_alignment = order_children_tips | set([x for i,x in outgroups.iteritems()])
                
                alignment_pruned = raxml_alignment.clone(2)
                
                for taxon in alignment_pruned.taxon_namespace:
                    if taxon.label not in tips_in_alignment:
                        alignment_pruned.remove_sequences([taxon])
                alignment_pruned.purge_taxon_namespace()
                alignment_pruned = remove_gaps(alignment_pruned)
                    
                #tax statements and constraints for families
                #constrain the current order
                for out_order, seq in outgroups.iteritems():
                    if out_order != order:
                        tax_statements.append('taxset ' + out_order + ' = ' + seq)
                    else:
                        tax_statements.append('taxset ' + order + ' = ' + ' '.join(order_children_tips))
                        con_statements.append('constraint ' + order + ' hard = ' + order)
                        
                #family-level constraints only if families represented in rainford tree        
                if order_children_families:  
                    #constrain monophyletic families
                    for family in set([x['family'] for i,x in align_data.iterrows() if x['raxml_placement'] == order and x['family'] is not np.nan]):
                        family_tips = set([x['genus'] + '_' + str(x['g_ott_id']) for i,x in align_data.iterrows() if x['family'] == family])
                        tax_statements.append('taxset ' + family + ' = ' + ' '.join(family_tips))
                        if len(family_tips) > 1:
                            if family in rainford_reverse.keys() or update_taxon(family) is not None:#if family in Rainford dictionary, it was already tested as monophyletic
                                con_statements.append('constraint ' + family + ' hard = ' +  ' '.join(family_tips))
                
                    #partial constraints on inter-familial relationships
                    constraint_idx = 1
                    for node in tree.postorder_internal_node_iter():
                        #positive part of constraint
                        children_families = set([leaf.taxon.label for leaf in node.leaf_nodes()])
                        
                        #negative part of constraint
                        families_not = [leaf.taxon.label for leaf in tree.leaf_nodes() if leaf.taxon.label not in children_families]
                        
                        if families_not:
                            con_statement = 'constraint ' + 'cons_' + str(constraint_idx) + ' partial = ' + \
                                            ' '.join(sorted(children_families)) + \
                                            ' : ' + \
                                            ' '.join(sorted(families_not))
                            constraint_idx += 1
                            con_statements.append(con_statement)
                
                #calibrations and constraints between orders
                constraint_idx = 1
                for node in order_tree.preorder_internal_node_iter():
                    children_orders = [leaf.taxon.label for leaf in node.leaf_nodes()]
                    
                    if set(children_orders) == set(['Ephemeroptera', 'Odonata']) and palaeoptera_paraphyletic is True:
                        #if palaeoptera is to be treated as paraphyletic, only add constraints and ignore calibrations
                        orders_external = ['Collembola', 'Protura', 'Diplura', 'Zygentoma', 'Archaeognatha']                    
                        
                        #constraint for Odonata + Neoptera
                        con_statement = 'constraint ' + 'order_cons_' + str(constraint_idx) + ' hard = ' + ' '.join(sorted(set(subordinal_trees.keys()) - set(orders_external) - set(['Ephemeroptera'])))
                        constraint_idx += 1
                        #constraint for Ephemeroptera + Odonata + Neoptera
                        con_statement = 'constraint ' + 'order_cons_' + str(constraint_idx) + ' hard = ' + ' '.join(sorted(set(subordinal_trees.keys()) - set(orders_external)))
                        constraint_idx += 1
                        continue
                        
                        
                    else:
                        con_statement = 'constraint ' + 'order_cons_' + str(constraint_idx) + ' hard = ' + ' '.join(sorted(children_orders))
                    
                    con_statements.append(con_statement)
                    
                    if order_tree.annotations.get_value('source') == 'rainford':
                        #to get calibration statements, calibrations have to be searched in the full rainford tree
                        children_tips = set()
                        for i_order in children_orders:
                            if subordinal_trees[i_order] is None:
                                children_tips.add(i_order)
                            else:
                                children_tips = children_tips | set([leaf.taxon.label for leaf in subordinal_trees[i_order].leaf_nodes()])
                                
                        rainford_names = [rainford_reverse[family] for family in children_tips]
                        node_rainford = rainford_tree.mrca(taxon_labels=rainford_names)
                        cal_mean = float(node_rainford.annotations.values_as_dict()['age_median'])
                    
                        
                    elif order_tree.annotations.get_value('source') == 'misof':
                        #for the misof tree, calibrations are already there as annotations
                        try:
                            cal_mean = float(node.annotations.get_value('age_median'))
                        except AttributeError:
                            cal_mean = []
                        
                    if cal_mean:    
                        cal_statement = 'calibrate ' + 'order_cons_' + str(constraint_idx) + ' = ' + \
                                'fixed(' + '{:.2f}'.format(cal_mean) + ')'
                        cal_statements.append(cal_statement)
                    else:
                        print 'No calibration for the ancestor of: ' + ','.join(children_orders)
                        
                    constraint_idx += 1 
                
                    
                #now write mrbayes file
                nexus_string = alignment_pruned.as_string(schema = 'nexus', unquoted_underscores=True)
                nexus_splits = nexus_string.split('BEGIN SETS;')
                mrbayes_out.write(nexus_splits[0])
                
                #now, start mrbayes block and add character partitions with GTR + GAMMA + I model
                mrbayes_out.write('\n\n' + 'BEGIN MRBAYES;' +  '\n\n')
                mrbayes_out.write(nexus_splits[1].split('END')[0] + '\n')
                mrbayes_out.write('partition by_gene = ' + 
                                str(len(alignment_pruned.character_subsets)) + 
                                ':' + 
                                ', '.join(alignment_pruned.character_subsets.keys()) +
                                ';\n')
                mrbayes_out.write('set partition = by_gene;\n')
                mrbayes_out.write('lset applyto=(all) nst=6 rates=invgamma;\n')
                mrbayes_out.write('unlink statefreq=(all) revmat=(all) shape=(all) pinvar=(all);\n')
                mrbayes_out.write('prset ratepr = variable;\n\n')
        
                #taxset, constraints and calibrations     
                #add tree and calibration priors
                constraints = [x.split(' ')[1] for x in con_statements]
                mrbayes_out.write('prset brlenspr = clock:birthdeath;\n')
                mrbayes_out.write(';\n'.join(tax_statements) + ';\n')
                mrbayes_out.write('\n\n')
                mrbayes_out.write(';\n'.join(con_statements)+ ';\n')
                mrbayes_out.write('\n\n')
                mrbayes_out.write(';\n'.join(cal_statements)+ ';\n')
                mrbayes_out.write('\n\n')                  
                mrbayes_out.write('prset topologypr = constraints(' + ','.join(constraints) + ');\n') 
                mrbayes_out.write('link topology=(all) brlens=(all);\n')
                
                
                priors = list()
                
                priors.append('sampleprob = ' + '{:.6f}'.format(len(alignment_pruned.taxon_namespace)/1000000.0)) #we sampled about N insect species form about 1M
                priors.append('speciationpr = exp(13.5)') #from the model with no shifts in condamine et al (2016) scientific reports (https://www.nature.com/articles/srep19208)
                priors.append('extinction = beta(2,2000)') #this puts extinction rate on the same order of magnitude as inferred in their models
                priors.append('nodeagepr = calibrated')
                priors.append('clockratepr = exponential(1000)') 
                priors.append('clockvarpr = igr') 
                priors.append('igrvarpr = exp(50)') 
                
                
                
                mrbayes_out.write('prset ' + ' '.join(priors) + ';\n')
                mrbayes_out.write('\n\n')
                        
                #finally, write mcmc command
                mrbayes_out.write('showmodel;\n')
                mrbayes_out.write('mcmcp ngen=500000000 diagnfreq=10000 relburnin=yes burninfrac=0.25 printfreq=50000 samplefreq=10000 checkfreq=100000 savebrlens=yes;\n')
                mrbayes_out.write('[mcmcp append=yes;]\n') #uncomment later in the mrbayes code if restart is needed
                mrbayes_out.write('mcmcp filename = ' + order + ';\n')
                mrbayes_out.write('mcmc; sumt; sump; \n\n END;')
    
            

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-a','--alignment', help = 'concatenated alignment in phylip format', default='ML_tree/genus_concatenated.phy')
    parser.add_argument('-p','--partitions', help = 'raxml partition file corresponding to the alignemnt', default='ML_tree/genus_partitions')
    parser.add_argument('-d','--aligndata', help = 'table with data about terminals in alignment', default='dataset_genus_overlap.csv')
    parser.add_argument('-f','--famtree', help = 'tree with family level constraints (rainford et al 2012)', default='constraint_info/rainford_family_tree.nxs')
    parser.add_argument('-i','--faminfo', help = 'table with taxon translation for rainford\'s family tree', default='constraint_info/rainford_family_key.txt')
    parser.add_argument('-o','--ordertree', help = 'tree with order level constraints (misof et al)', default='constraint_info/fully_annotated_misof.nexml')
    
  
    args = parser.parse_args()
    #args = parser.parse_args([])
    
    print 'Opening input files'
    #opening input files and doing some basic data preparation
    rainford_trees = dendropy.TreeList.get(path = args.famtree,schema="nexus", preserve_underscores=True)
    rainford_tree = rainford_trees[1] #this gets the dated tree
    rainford_tree.annotations.add_new('source', 'rainford') 
    
    misof_tree = dendropy.Tree.get(path = args.ordertree, schema = 'nexml')
    misof_tree.annotations.add_new('source', 'misof')


    align_data = pandas.read_csv(args.aligndata, sep='\t')
    
    raxml_alignment = dendropy.DnaCharacterMatrix.get(path = args.alignment, schema = 'phylip')
    
    raxml_partitions = open(args.partitions,'r')
    for line in raxml_partitions:
        name = line.split('=')[0].split(' ')[-1]
        idx_range = line.split('=')[-1].split('-')
        raxml_alignment.new_character_subset(label = name, character_indices = range(int(idx_range[0]) - 1,int(idx_range[1])))
    raxml_partitions.close()
    
    rainford_df = pandas.read_csv(args.faminfo, sep ="\t", header=None)
    rainford_df.columns = ['rainford_order','rainford_family','rainford_tipname']   
    #update order names 
    
    
    #make several dicts to relate rainford tip names, rainford families and their currently accepted names in OTL (dropping non-monophyletic families)
    
    rainford_key = {x[2]:x[1] for i,x in rainford_df.iterrows()}
    rainford_taxa = [x[1] for i,x in rainford_df.iterrows()]
    print 'Checking OTL for monophyly of tips in Rainford\'s tree'
    rainford_updated = [update_taxon(x) for x in rainford_taxa]
    rainford_df['otl_taxon'] = rainford_updated 
    rainford_key_updated = {x[2]:x['otl_taxon'] for i,x in rainford_df.iterrows()}
    rainford_reverse = {x:i for i,x in rainford_key_updated.iteritems() if x is not None}
    #rainford_df.to_csv('rainford_family_otl.txt', sep = '\t')
    
    print 'Pruning Rainford\'s tree to families (monophyletic on OTL) present in the egg dataset and sequence data'
    #prune the rainford tree to families in our egg dataset and alignment
    all_egg_taxa = set(align_data.family).union(align_data.raxml_placement)
    egg_in_constraint = [rainford_reverse[x] for x in all_egg_taxa if x in rainford_reverse.keys()]
    pruned_tree = rainford_tree.clone(2).extract_tree_with_taxa_labels(egg_in_constraint)
    pruned_tree.update_bipartitions(suppress_unifurcations=False)
    #rename taxa in pruned_tree
    for leaf in pruned_tree.leaf_node_iter():
        leaf.taxon.label = rainford_key_updated[leaf.taxon.label]

    #now check if all orders still represented
    tip2order = {x.family:x.raxml_placement for i,x in align_data.iterrows()}
    tip2order.update({x.raxml_placement:x.raxml_placement for i,x in align_data.iterrows()})
    orders_in_tree = sorted(set([tip2order[leaf.taxon.label] for leaf in pruned_tree.leaf_node_iter()]))
    orders_not_represented = set(align_data.raxml_placement) - set(orders_in_tree)
    
    if orders_not_represented:
        print 'The following orders have no monophyletic families in our dataset, only order-level constraint used:'
        print ', '.join(sorted(orders_not_represented))
    #if not all orders represented, first collapse the tree to keep this order
    if orders_not_represented:
        for order in orders_not_represented:
            pruned_tree = rainford_tree.clone(2)
            rainfordtips_for_this_order = set(rainford_df.loc[rainford_df.rainford_order == order,'rainford_tipname'])
            mrca_node = rainford_tree.mrca(taxon_labels = rainfordtips_for_this_order)
            mrca_tips = mrca_node.leaf_nodes()
            keep = mrca_tips[0]
            keep.taxon.label = order
            pruned_tree.update_taxon_namespace()
            egg_in_constraint.append(order)
            rainford_reverse[order] = order
            rainford_key_updated[order] = order
            
        #Now, prune to families/orders in the dataset again
        pruned_tree = pruned_tree.extract_tree_with_taxa_labels(egg_in_constraint).clone(2)
        pruned_tree.update_bipartitions(suppress_unifurcations=False)
        #rename taxa in pruned_tree
        for leaf in pruned_tree.leaf_node_iter():
            leaf.taxon.label = rainford_key_updated[leaf.taxon.label]
        
    
    #separate an order level tree and subordinal trees
    #print the order-level tree to make sure everything is fine
    #rainford used order names when they only had one family per order
    print 'Separating tree in backbone order-level and within-order trees'
    rainford_order_tree = pruned_tree.clone(2)
    subordinal_trees = dict()
    for order in set(align_data.raxml_placement):
        tip_labels = [leaf.taxon.label for leaf in rainford_order_tree.leaf_node_iter()]
        if order in tip_labels: #if order is a tip label, do nothing and only record that it exists
            subordinal_trees.update({order:None})
        else:    
            families = set(align_data.loc[[x == order for x in align_data.raxml_placement],'family'])
            families_in_tree = families.intersection(tip_labels)
            if families_in_tree:
                #to test monophyly, all descendants of the mrca have to be these families
                mrca_node = rainford_order_tree.mrca(taxon_labels=families_in_tree)
                mrca_descendants = set([leaf.taxon.label for leaf in mrca_node.leaf_iter()])
                if families_in_tree == mrca_descendants:
                    subordinal_trees.update({order:mrca_node.extract_subtree()})
                    mrca_node.clear_child_nodes()
                    mrca_node.taxon = dendropy.Taxon(order)
            
    rainford_order_tree.update_taxon_namespace()
    rainford_order_tree.purge_taxon_namespace()
    rainford_order_tree.update_bipartitions(suppress_unifurcations=False)  
    rainford_order_tree.annotations.add_new('source', 'rainford')
    print 'Rainford backbone tree used for constraints:'      
    rainford_order_tree.print_plot()
    rainford_order_tree.write(path='mrbayes_smalltrees/rainford_backbone.nexml', schema='nexml')
    rainford_order_tree.write(path='mrbayes_smalltrees/rainford_backbone.tre', schema='nexus')
    
    #now make a misof backbone:
    print 'Pruning Misof tree to one tip per order in the dataset'
    #since our misof tree already has nodes with order names, getting the backnone is easy
    #we only prune all but one children of nodes that have a property named 'order'
    tips_removed = True
    while tips_removed:
        tips_removed = False
        for node in misof_tree.postorder_node_iter():
            if 'order' in node.annotations.values_as_dict().keys() and \
               node.annotations.get_value('order') in set(align_data['raxml_placement']) and \
               (node.taxon is None or \
               node.taxon.label != node.annotations.get_value('order')):
                                   
                node.taxon = dendropy.Taxon(label = node.annotations.get_value('order'))
                if node.num_child_nodes():
                    node.clear_child_nodes()
                tips_removed = True
                break
    misof_tree.update_taxon_namespace()
                    
    tips_to_keep = [tip.taxon.label for tip in misof_tree.leaf_node_iter() if tip.taxon.label in set(align_data['raxml_placement'])]
    
    misof_tree.retain_taxa_with_labels(tips_to_keep, update_bipartitions=True, suppress_unifurcations=False)
    misof_tree.update_taxon_namespace()
    misof_tree.purge_taxon_namespace()

    print 'Misof backbone tree used for constraints:' 
    misof_tree.print_plot()
    
    
    #To make it easier when stitching, transform tree branch lengths to correspond to node ages
    for node in misof_tree.postorder_node_iter():
        if node.is_leaf():
            node.edge_length = float(node.parent_node.annotations.get_value('age_median'))
        elif node != misof_tree.seed_node:
            node.edge_length = float(node.parent_node.annotations.get_value('age_median')) - float(node.annotations.get_value('age_median'))
    misof_tree.calc_node_ages()       
        
        
        
    #stitching is failing later because there is a problem with tree rooting, this should save the backbone appropriately
    misof_tree.reroot_at_node(misof_tree.seed_node, update_bipartitions=False, suppress_unifurcations=False)

    #save tree
    misof_tree.write(path='mrbayes_smalltrees/misof_backbone.nexml', schema='nexml')
    misof_tree.write(path='mrbayes_smalltrees/misof_backbone.tre', schema='nexus')
    
    
    
    print 'Making mrbayes files with genera as tips'  
    #small trees with misof do not make much sense, since we do not have an ultrametric backbone tree
    #Misof et al do not provide the ultrametric tree, so we only have topology and the ages reported in their supplement
    #however, not all nodes had reported ages
    
    print "1 - Rainford family and order backbone, Rainford calibration"
    make_big_tree(raxml_alignment, align_data, rainford_order_tree, subordinal_trees, outfolder='mrbayes_bigtree/rainford', palaeoptera_paraphyletic = False)
    make_orders_trees(raxml_alignment, align_data, rainford_order_tree, subordinal_trees, outfolder='mrbayes_smalltrees/rainford', palaeoptera_paraphyletic = False)
    
    print "2 - Misof et al's order backbone and calibration, rainford family backbone"
    make_big_tree(raxml_alignment, align_data, misof_tree, subordinal_trees, outfolder='mrbayes_bigtree/misof', palaeoptera_paraphyletic = False)
    make_orders_trees(raxml_alignment, align_data, misof_tree, subordinal_trees, outfolder='mrbayes_smalltrees/misof', palaeoptera_paraphyletic = False)
    
    print "3 - Misof et al's order backbone and calibration, rainford family backbone, palaeoptera paraphyletic"
    make_big_tree(raxml_alignment, align_data, misof_tree, subordinal_trees, outfolder='mrbayes_bigtree/misof_palaeoptera_paraphyletic', palaeoptera_paraphyletic = True)
    #make_orders_trees(raxml_alignment, align_data, misof_tree, subordinal_trees, outfolder='mrbayes_smalltrees/misof_palaeoptera_paraphyletic', palaeoptera_paraphyletic = True)

    print 'DONE'
