#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 23 17:24:27 2019

@author: Bruno de Medeiros

This script calculates the standard deviation of split frequencies for all 
trees inferred. 
Since it takes a long time, we run it using ipyparallel, writing the most
time-consuming operations in a function

"""

import argparse, time, ipyparallel as ipp, sys
    
def get_stdev_split_freqs(order, backbone):
    import sys, dendropy, gzip, math, glob #os, psutil
    #p = psutil.Process(os.getpid())
    #p.cpu_affinity([cpu])
    
    sys.stderr.write('Working on ' + order + ' ' + backbone + '\n')
    sys.stderr.flush()

    def read_one_file(tree_file, taxa, p_burnin):
        sys.stderr.write('Reading: ' + tree_file + '\n')
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
        
        # calculate tre
        burnin = math.ceil(p_burnin * ntrees)
        
        # read trees except for burnin
        tr = dendropy.TreeArray(taxon_namespace = taxa, 
                                ignore_edge_lengths = True,
                                ignore_node_ages = True,
                                is_rooted_trees = True)
        tr.read(path = tree_file,
                              schema = 'nexus',
                              tree_offset=burnin)
        
        sys.stderr.write('Finished: ' + tree_file + '\n')
        return tr.split_distribution.calc_freqs()    


    # find tree pahts
    tree_paths = glob.glob('mrbayes_smalltrees/treelinks_' + backbone + '/' + order + '*.t')

    
    # read first tree to get a taxon namespace
    infile = open(tree_paths[0],'r')
    tree_st = []
    for i, line in enumerate(infile):
        tree_st.append(line)
        if 'tree gen' in line:
            tree_st.append(line)
            break
    infile.close()
    taxa = dendropy.Tree.get(data = '\n'.join(tree_st), schema='nexus').taxon_namespace
    

    split_freqs = [read_one_file(tree_path,taxa,0.1) for tree_path in tree_paths]

    
    #let's use pandas to calculate stdev of split frequencies
    import pandas as pd
    outdict = {'backbone': backbone,
            'order': order,
            'stdev': pd.DataFrame(split_freqs).fillna(0).std(0).mean(0)}
    sys.stdout.write(str(outdict) + '\n')
    sys.stdout.flush()
    return outdict

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-p','--profile',help='ipcluster profile to use') 
    parser.add_argument('-i','--id',help='ipcluster id to use') 
 
    args = parser.parse_args()
    
    
    # we will loop through backbones and orders
    backbones = ['misof', 'rainford'] 
    orders = ['Archaeognatha',
                'Coleoptera',
                'Collembola',
                'Dermaptera',
                'Dictyoptera',
                'Diplura',
                'Diptera',
                'Ephemeroptera',
                'Grylloblattodea',
                'Hemiptera',
                'Hymenoptera',
                'Lepidoptera',
                'Mantodea',
                'Mantophasmatodea',
                'Mecoptera',
                'Megaloptera',
                'Neuroptera',
                'Odonata',
                'Orthoptera',
                'Phasmatodea',
                'Plecoptera',
                'Protura',
                'Psocodea',
                'Siphonaptera',
                'Strepsiptera',
                'Thysanoptera',
                'Trichoptera',
                'Zygentoma']
    
    bb_expanded = list()
    or_expanded = list()
    
    for backbone in backbones:
        for order in orders:
            bb_expanded.append(backbone)
            or_expanded.append(order)
    

    #set up parallel run
    rc = ipp.Client(profile=args.profile, cluster_id=args.id)
    dview = rc[:]
    #dview.scatter("cpu", range(len(dview)), flatten=True)

    asdf = dview.map_async(get_stdev_split_freqs, 
                                 or_expanded, 
                                 bb_expanded)

    
    
    #print messages while not ready
    (stdout0, stderr0) = (['']*len(dview), ['']*len(dview))
    while not asdf.ready():
    # check if stdout changed for any kernel
        if (asdf.stdout, asdf.stderr) != (stdout0, stderr0):
            for i in range(len(asdf.stderr)):
                if asdf.stderr[i] != stderr0[i]: 
                    # print only new stdout's without previous message
                    previous = len(stderr0[i])
                    stderr0[i] = asdf.stderr[i]
                    sys.stderr.write('kernel ' + str(i) + ': ' + stderr0[i][previous:])
                elif asdf.stdout[i] != stdout0[i]:
                    previous = len(stdout0[i])
                    stdout0[i] = asdf.stdout[i]
                    sys.stdout.write('kernel ' + str(i) + ': ' + stdout0[i][previous:])
            sys.stderr.flush()
            sys.stdout.flush()
        else:
            time.sleep(10)
        
    
    
    
            
    #save as csv
    import pandas as pd
    pd.DataFrame(asdf.get()).to_csv('final_trees/stdev_split_freqs.csv')
            
    
    
