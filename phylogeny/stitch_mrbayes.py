#!/usr/bin/env python

#This script gets mrbayes results for each order and grafts them onto the backbone tree
#This is done both for the MCC tree and for a sample of the posterior distribution


import dendropy, pandas, argparse, glob, random,gzip,multiprocessing
from math import ceil
from os.path import basename
from joblib import Parallel, delayed



def brlen2time(tree):
    clockrate = float(tree.annotations.get_value('clockrate'))
    for node in tree.postorder_node_iter():
        if node.edge_length:
            node.edge_length = node.edge_length / clockrate
    
    return tree
    

#This function finds the tree indexes that have to be read:
def count_trees(tree_path, burnin, n_per_treefile):
    tree_file = tree_path
    ntrees = 0
    if '.gz' in tree_file:
        infile = gzip.open(tree_file, 'rb')
    else:
        infile = open(tree_file,'r')
        
    first_tree = 0
    for i, line in enumerate(infile):
        if 'tree gen' in line:
            ntrees += 1
            if not first_tree: first_tree = i
        
    infile.close()
                
    n_burnin = int(ceil(float(ntrees) * burnin))
    step_size = max([len(range(n_burnin+1,ntrees))/n_per_treefile,1])
    tree_idx = range(n_burnin,ntrees,step_size)[-n_per_treefile:]
    
        
    return tree_idx


#This function reads trees based on predefined indexes and a taxon namespace
def read_trees(tree_path,nfiles,taxa):
    
    tree_indices = count_trees(tree_path, burnin=0.1, n_per_treefile=10000/nfiles)
    
    print 'Reading ' + str(len(tree_indices)) + ' posterior trees from ' + tree_path + '.'
    
    #to speed up tree reading, we will first work with strings and then read trees
    
    tree_strings = []
    if '.gz' in tree_path:
        infile = gzip.open(tree_path,'rb') 
    else:
        infile = open(tree_path,'r') 
        

    lines_pre_trees = 0   
    for i,line in enumerate(infile):
        if 'tree gen' not in line:
            tree_strings.append(line)
            lines_pre_trees += 1
            
        elif i - lines_pre_trees in tree_indices:
                tree_strings.append(line)
                            
    tree_strings.append('END;')
    
    infile.close()
            
        
    
    trees = dendropy.TreeList.get(data=''.join(tree_strings), schema='nexus', preserve_underscores=True, taxon_namespace = taxa)
    return trees

            


#function updates taxon names in rainford according to open tree of life, and returns None for non-monophyletic groups
#monophyly is assumed if a node exists in open tree of life
def get_mcc_and_posterior_order(tree_paths, focal_order, taxdata, burnin=0.1, nmax = 10000):
    #let's first create the taxon namespace
    if '.gz' in tree_paths[0]:
        handle = gzip.open(tree_paths[0],'rb')
    else:
        handle = open(tree_paths[0],'rb')
    
    temp_tree = dendropy.Tree.get(file = handle, schema='nexus', preserve_underscores=True)
    handle.close()
    taxa = temp_tree.taxon_namespace
        
    
    
    #Now that we have the tree indexes and a taxon name space, let's load all trees
    
    
    #I tried to parallelize this part, but joblib cannot take nested parallel loops
    #all_post_trees = Parallel(n_jobs=min(threads,len(tree_paths)),backend = 'threading', verbose=50) \
    #                         (delayed(read_trees)(tree_path,tree_idx,taxa) for tree_path in tree_paths)
    #Using a regular loop instead
    
    all_post_trees = [read_trees(tree_path,len(tree_paths),taxa) for tree_path in tree_paths]
    
    
    post_trees = dendropy.TreeList()
    
    for posterior in all_post_trees:
        post_trees.extend(posterior)
    
    #now, rescale branch lengths to units of time
    for tree in post_trees:
        tree = brlen2time(tree)
        
    #now get mcc tree 
    
    print 'Finding MCC tree for ' + focal_order
    #something wrong with sumtrees.py, forget about branch lengths for now
    #post_trees.write(path='.temp_tree', schema='nexus')
    #sumtree_command = 'sumtrees.py --ultrametric --summary-target mcct --burnin 0 --edges median-age -M -F nexus -i nexus -v 0 .temp_tree'
    #mcc_tree = dendropy.Tree.get(data=subprocess.check_output(sumtree_command.split()), schema='nexus')
   
    mcc_tree = post_trees.maximum_product_of_split_support_tree()
    
    
    print 'Extracting ingroup for ' + focal_order
    #now prune to this order only
    order_children_tips = set([x['genus'] + '_' + str(x['g_ott_id']) for i,x in taxdata.iterrows() if x['raxml_placement'] == focal_order])
    
    
    tips_to_keep = set([leaf.taxon for leaf in mcc_tree.leaf_node_iter() \
                     if leaf.taxon.label in order_children_tips])
                     
    mcc_mrca_node = mcc_tree.mrca(taxa=tips_to_keep)
    mcc_mrca_node.annotations.add_new('misof_clade',focal_order)
    
    #and do the same for all posterior trees:
    post_mrca_nodes = []
    for tree in post_trees:
        temp_node = tree.mrca(taxa = tips_to_keep)
        temp_node.annotations.add_new('misof_clade',focal_order)
        post_mrca_nodes.append(temp_node)
    
    return {'mcc':mcc_mrca_node,'posterior':post_mrca_nodes}

def get_subtree(order, small_trees, align_data, backbone_tree):
    print 'Retrieving subtree for ' + order
    
    backbone_tree = backbone_tree.clone(2)
    
    tips_in_order = [x['genus'] + '_' + str(x['g_ott_id']) for i,x in align_data.iterrows() if x['raxml_placement'] == order]
    
    #at this point, there shouldn't be any order in the tree without species in the alignment
    if not tips_in_order:
        raise Exception('order in backbone but not in alignment')
        
    #if only one tip for an order, simply rename it    
    elif len(tips_in_order) == 1:
        print order + ' has only one tip. Substituting genus for order.'
        order_node = backbone_tree.find_node_with_taxon_label(order)
        order_node.taxon.label = tips_in_order[0]
        order_node.annotations.add_new('misof_clade',order)
        
        mcc_and_posterior = {'mcc':order_node,'posterior':[order_node]}            
    
    #if more than one tip, replace with mrbayes tree
    else:
        mcc_and_posterior = get_mcc_and_posterior_order(tree_paths = small_trees[order],
                                                        focal_order = order,
                                                        taxdata = align_data,
                                                        burnin = 0.1,
                                                        nmax = 10000)
                                                        
    return {order:mcc_and_posterior}



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-d','--aligndata', help = 'table with data about terminals in alignment', default='dataset_genus_overlap.csv')
    parser.add_argument('-b','--backbone', help = 'backbone order-level tree in nexml format', default='mrbayes_smalltrees/rainford_backbone.nexml')
    parser.add_argument('-n','--n-posterior',help = 'number of posterior trees to produce', default = 100, type=int)
    parser.add_argument('-T','--threads',type=int, help = 'number of threads, defaults to available cores')
    parser.add_argument('-P','--parallel',action = 'store_true', help = 'whether to parse subtrees in parallel (faster, but consumes a lot of memory)')
    parser.add_argument('-t','--trees',help = 'list of order-level trees obtained from mrbayes', default = 'mrbayes_smalltrees/*.t.gz', nargs='+')
    
    
    args = parser.parse_args()
    #args = parser.parse_args('--threads 1 --aligndata ../../dataset_genus_overlap.csv --n-posterior 5 --backbone ../rainford_backbone.nexml --trees *.t'.split())
    
    if not args.threads:
        args.threads = multiprocessing.cpu_count()
    
    print 'Opening input files'
    #opening input files and doing some basic data preparation
    order_tree = dendropy.Tree.get(path = args.backbone,
                                   schema="nexml")
    
    align_data = pandas.read_csv(args.aligndata, sep='\t')
    
    small_trees = dict()
    
    if '*' in args.trees[0]: #if running default, has to glob
        args.trees = glob.glob(args.trees[0])
    else: #otherwise, there should be already a list of files
        pass
        
    #register which small trees correspond to which order    
    for small_tree_path in args.trees:
        order = basename(small_tree_path).split('.')[0]
        if order in small_trees.keys():
            small_trees[order].append(small_tree_path)
        else:    
            small_trees[order] = [small_tree_path]
    
    print 'Stitching small trees onto backbone'
    
    #first, clone order_tree to start a final mcc tree and N posterior trees
    final_mcc_tree = order_tree.clone(2)
    final_posterior = dendropy.TreeList([order_tree.clone(2) for i in xrange(args.n_posterior)])
    
    #get all orders in backbone tree
    orders = [leaf.taxon.label for leaf in final_mcc_tree.leaf_node_iter()]
    
    #get subtrees for all orders
    if args.parallel:
        subtrees = Parallel(n_jobs=min(args.threads,len(orders)), verbose=50) \
                             (delayed(get_subtree)(order, small_trees, align_data, final_mcc_tree) for order in orders)
        subtree_dict = dict()
        for sub in subtrees:
            subtree_dict.update(sub)
    
    #for each order, replace tip for the mcc tree for the order
    for order in orders:
        if args.parallel:
            order_trees = subtree_dict[order]
        else:
            order_trees = get_subtree(order, small_trees, align_data, order_tree)[order]
            
        print 'Stitching subtrees for ' + order

        #MCC
        print 'Stitching MCC'                                                
        order_leaf = final_mcc_tree.find_node_with_taxon(lambda x: x.label == order)
        order_parent = order_leaf.parent_node
        
        order_parent.add_child(order_trees['mcc'])
        order_parent.remove_child(order_leaf, suppress_unifurcations = False)
        
        #POSTERIOR
        print 'Stitching ' + str(args.n_posterior) + ' posterior trees'
        for tree in final_posterior:
            order_leaf = tree.find_node_with_taxon(lambda x: x.label == order)
            order_parent = order_leaf.parent_node

            order_parent.add_child(random.choice(order_trees['posterior']))
            order_parent.remove_child(order_leaf, suppress_unifurcations = False)
            
        

        final_mcc_tree.update_taxon_namespace()
        final_mcc_tree.update_bipartitions(suppress_unifurcations=False)
        final_posterior.update_taxon_namespace()
        for tree in final_posterior:
            tree.update_bipartitions(suppress_unifurcations=False)
        
    
    
    final_mcc_tree.write(path = 'mcc_tree.tre', schema = 'nexus', unquoted_underscores=True)
    final_posterior.write(path = 'posterior_trees.tre', schema = 'nexus', unquoted_underscores=True)
    print 'DONE'       
