#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 29 2022

@author: Palash Sashittal
"""

import pandas as pd
import sys
import argparse
import itertools
import math
import numpy as np
import re
from scipy.stats import betabinom

import networkx as nx

def writeDOT(T, dot_file):
    with open(dot_file, 'w') as output:

        output.write(f'digraph N {{\n')
        output.write(f"\toverlap=\"false\"\n")
        output.write(f"\trankdir=\"TB\"\n")

        idx_dict = {}
        idx = 0
        for node in T.nodes:
            idx_dict[node] = idx
            output.write(f'\t{idx} [label=\"{node}\", style=\"bold\"];\n')
            idx += 1
        
        for edge in T.edges:
            output.write(f"\t{idx_dict[edge[0]]} -> {idx_dict[edge[1]]} [style=\"bold\"];\n")
        
        output.write(f'}}')    

def get_loc(header):
  pattern = r"chr(\d+|X|Y)-(\d+)"
  match = re.search(pattern, header)
  if match:
    chr = match.group(1)
    offset = match.group(2)
    return (header, chr, offset)
  else:
    return None

def write_bed_file(df, output_file):
  cols = df.columns
  col_list = list(cols)
  col_list.remove('cluster_id')
  headers = [s[:-2] for s in col_list]

  loc_data = []

  for elt in headers:
    tup = get_loc(elt)
    loc_data.append(tup)

  sorted_loc_data = sorted(loc_data, key=lambda x: (int(x[1]), int(x[2])))

  sorted_columns = [tup[0] for tup in sorted_loc_data]

  with open(output_file, "w") as file:
      for tup in sorted_loc_data:
          file.write("{chr}\t{offset}\t{offset}\t{name}\n".format(name=tup[0], chr=tup[1], offset=tup[2]))


def write_nexus_file(df, output_file):
  df1 = df.drop('cluster_id', axis=1)
  cols = df1.columns
  labels = list(cols)
  cells = list(df1.index)
  cells.append('O')

  with open(output_file, "w") as file:
      file.write("#NEXUS\n\n")
      file.write("begin data;\n")
      file.write("\tDimensions: ntax={cols}, nchar={rows};\n".format(cols=len(cells), rows=len(labels),))
      file.write("\tFormat datatype=standard missing=? gap=-;\n")
      file.write("\tMatrix\n")

      for cell in cells:
          file.write("{cell} ".format(cell=cell))
          for label in labels:
              if(cell == 'O'):
                 value = 0
                 file.write("{value}".format(value=value))
              else:
                 value = df1.loc[cell, label]
                 file.write("{value}".format(value=value))

          file.write("\n")

      file.write("\t;\n")
      file.write("End;\n")


def tree_to_newick(T, root=None):
    if root is None:
        roots = list(filter(lambda p: p[1] == 0, T.in_degree()))
        assert 1 == len(roots)
        root = roots[0][0]
    subgs = []
    while len(T[root]) == 1:
        root = list(T[root])[0]
    for child in T[root]:
        while len(T[child]) == 1:
            child = list(T[child])[0]
        if len(T[child]) > 0:
            child_newick = tree_to_newick(T, root=child)
            if child_newick != '()':
                subgs.append(child_newick)
        else:
            if child.startswith('s'):
                subgs.append(child)
    # return "(" + ','.join(map(str, subgs)) + ")"
    if len(subgs) == 1:
        return str(subgs[0])
    else:
        return "(" + ','.join(map(str, subgs)) + ")"    

def dfs_with_path(graph, start, target, path=None, visited=None):
    if path is None:
        path = []
    if visited is None:
        visited = set()

    path.append(start)
    visited.add(start)

    if start == target:
        return path

    for neighbor in graph.successors(start):
        if neighbor not in visited: 
            result = dfs_with_path(graph, neighbor, target, path.copy(), visited)
            if result: 
                return result
    return None

def main(args):

 
    print("HELLO WORLD!!!!!!!!!!!!!!!!!!!!!")
    
    np.random.seed(args.s)

    T = nx.DiGraph() # mutation tree
    Tc = nx.DiGraph() # copy number tree

    # add root nodes
    T.add_node('root')
    Tc.add_node(0)

    # build tree
    ncharacters = args.m
    nmutations = ncharacters
    nclusters = args.p
    max_losses = args.k
    max_cn = args.maxcn
    mutation_rate = args.l
    nnodes = ncharacters + nclusters
    print(max_cn)

    character_list = [f'c{character_index}' for character_index in range(ncharacters)]
    cluster_list = [f'd{cluster_index}' for cluster_index in range(1, nclusters)]
    #event randomization
    event_order = np.random.permutation(character_list + cluster_list)
    #vector tracking losses per character
    loss_counter = np.zeros((ncharacters, 1))
    #dict. tracking which cluster index losses occur on
    loss_dictionary = {f'd{cluster_index}': [] for cluster_index in range(1, nclusters)}
    
    #one-to-one correspondence b/w events and characters
    B = np.zeros((ncharacters + nclusters, ncharacters + 1), dtype=int)
    #one-to-one correspondence b/w cn and characters
    R = np.zeros((nclusters, ncharacters), dtype=int)    
    #random generation of starting copy number
    #need to normalize this value
    R[0, :] = np.random.randint(max_cn - max_losses - 1, size = ncharacters) + max_losses + 1
    
    #for each event
    for node_index, event in enumerate(event_order):
        #how many mutations have occurred up to (not including) node_index
        nprev_mutations = sum([1 for x in event_order[:node_index] if x.startswith('c')])
        #one-indexing
        node_index += 1

        #if its a change in cluster id
        if event.startswith('d') and nprev_mutations > 0:
            #find a valid location for cluster change (parent must be a mutation)
            while parent_node.startswith('d') or parent_node == 'root':
                #constant attachment kernel for growing random network
                parent_node_index = np.random.randint(node_index)
                parent_node = list(T.nodes)[parent_node_index]
        #if its a mutation or there are no previous mutations
        else:
            #constant attachment kernel for growing random network
            parent_node_index = np.random.randint(node_index)
            parent_node = list(T.nodes)[parent_node_index]
            

        #add the random node/event to the network
        T.add_edge(parent_node, event)
        #add the random node/event to the matrix B, init equal to its parent
        B[node_index, :] = B[parent_node_index, :]
        #if its a change in cluster id
        if event.startswith('d'):
            #get the new cluster id
            cluster_id = int(event.lstrip('d'))
            #change the cluster id to the new one
            B[node_index, -1] = cluster_id
            #get the parent cluster id
            parent_cluster_id = B[parent_node_index, -1]
            #add edge from old to new cluster id in the copy number tree (for scarlet)
            Tc.add_edge(parent_cluster_id, cluster_id)
            #init the new copy number profile to parent values
            R[cluster_id, :] = R[parent_cluster_id, :]

            #Rather than deleting each mutation at random
            #new algo:
            #x = roll dice with k sides where k = number of mutations b/w node_index and root
            #delete the last x mutations in order
            #this makes sense because we will eventually impose a ordering to these mutations

            for mutation in range(ncharacters):
                #do dfs to get the mutations on path from root to mutation
                full_path = dfs_with_path(T,'root',event)
                full_path.reverse()

                #print(f"full path for mutation losses: {full_path}")

                x =np.random.randint(1,len(full_path))
                subpath = []
                


                for i in range(x+1):
                    subpath.append(full_path[i])
                


                print(f"subpath for mutation losses: {subpath}")
                
                print(B[parent_node_index])

                mut = "c"+str(mutation)

                if B[parent_node_index, mutation] == 1 and loss_counter[mutation] < max_losses: 
                    print("HERE")
                    #add mutation loss to event matrix
                    B[node_index, mutation] = loss_counter[mutation] + 2
                    #increment loss counter for this mutation
                    loss_counter[mutation] += 1
                    #add the mutation to the loss dictionary
                    loss_dictionary[event].append(mutation)
                    #decrement the copy number by 1 at this location
                    R[cluster_id, mutation] -= 1
        
        #if its a mutation gain
        elif event.startswith('c'):
            #add the mutation to the event matrix
            mutation = int(event.lstrip('c'))
            B[node_index, mutation] = 1
            #perhaps add additional check "is this a descendant of a CNV"
            #this is because the loss will be transferred to all children

        if args.v:
            print(parent_node_index, parent_node, node_index, event)
    

    # randomize the copy number states for mutations that have never been lost
    for mutation in range(nmutations):
        if loss_counter[mutation] == 0:
            for cluster_id in range(nclusters):
                R[cluster_id, mutation] = np.random.randint(max_cn - 1) + 1

    # check that all copy number states are non-zer positive
    assert(len(np.where(R == 0)[0]) == 0)
    
    # check all SNV losses are supported by CNVs
    for cn_edge in Tc.edges:
        for mutation in loss_dictionary[f'd{cn_edge[1]}']:
            assert(R[cn_edge[0], mutation] > R[cn_edge[1], mutation])    


    #simulate location data

    cnv_size = 5000000
    num_chr = 22
    chromosome_sizes = [
        247249719, 242951149, 199501827, 191273063, 180857866, 170899992,
        158821424, 146274826, 140273252, 135374737, 134452384, 132349534,
        114142980, 106368585, 100338915, 88827254, 78774742, 76117153,
        63811651, 62435964, 46944323, 49691432
    ]

    losses = []

    for i in range(nmutations):
      if(loss_counter[i] > 0):
        losses.append(i)

    chr = np.random.randint(1, num_chr+1)
    offset = np.random.randint((cnv_size/2), chromosome_sizes[chr-1]-(cnv_size/2))

    names = {}

    for i in range(len(losses)):
      name = f"c{losses[i]}-chr{chr}-{offset}"
      names[losses[i]] = name

    for i in range(1,nmutations+1):
      if(names.get(i) == None):
        randchr = np.random.randint(1, num_chr+1);
        randoffset = np.random.randint((cnv_size/2), chromosome_sizes[randchr-1]-(cnv_size/2))
        name = f"c{i}-chr{randchr}-{randoffset}"
        names[i] = name


    if args.v:
        print('-'*50)
        print('loss counter')
        print('-'*50)
        print(loss_counter)
        print('-'*50)
        print('loss dictionary')
        print('-'*50)
        print(loss_dictionary)
        print('-'*50)
        print('copy number states')
        print('-'*50)
        print(R)
            
    # assign cells and generate character-state matrix
    leaf_indices = []
    for idx, node in enumerate(T.nodes):
        if len(T[node]) == 0:
            leaf_indices.append(idx)    
    nleaves = len(leaf_indices)
    
    ncells = args.n
    assert(ncells > nclusters)
    cell_assignment = np.random.randint(ncharacters, size=ncells-nleaves)
    complete_cell_assignment = list(cell_assignment) + leaf_indices
    Bcell = B[complete_cell_assignment, :]
        
    # observed matrix
    A = B.copy()
    for mutation in range(ncharacters):
        A[A[:,mutation] > 1, mutation] = 0
    Acell = A[complete_cell_assignment, :]
    
    # cell tree
    celltree = T.copy()
    for cell_id, assigned_node_index in enumerate(complete_cell_assignment):
        celltree.add_edge(list(T.nodes)[assigned_node_index], f's{cell_id}')

        
    # generate read counts
    mean_coverage = args.cov
    fp_rate = args.a
    fn_rate = args.b
    ado_precision = args.ado

    Rtotal = np.zeros((ncells, nmutations), dtype=int)
    Vcount = np.zeros((ncells, nmutations), dtype=int)
    for cell in range(ncells):
        for mutation in range(nmutations):
            cluster_id = Acell[cell, -1]
            nvariant = Acell[cell, mutation]
            ntotal = R[cluster_id, mutation]

            latent_vaf = nvariant / ntotal

            nreads = np.random.poisson(mean_coverage)
            Rtotal[cell, mutation] = int(nreads)

            post_error_vaf = fp_rate + (1 - fp_rate - fn_rate) * latent_vaf
            ado_alpha = post_error_vaf * ado_precision
            ado_beta = ado_precision * (1 - post_error_vaf)
            nvariant_reads = betabinom.rvs(nreads, ado_alpha, ado_beta)

            Vcount[cell, mutation] = int(nvariant_reads)
            
    # generate the binarized mutation matrix
    vaf_threshold = args.vafthreshold
    variant_read_threshold = args.readthreshold
    VAF_mat = Vcount / Rtotal
    mutation_mat = ((VAF_mat >= vaf_threshold) & (Vcount >= variant_read_threshold)).astype(int)
    mutation_mat = np.hstack((mutation_mat, Acell[:,-1][:,np.newaxis]))
    
    # introduce missing entries
    Acell_missing = Acell.copy()
    Rtotal_missing = Rtotal.copy()
    Vcount_missing = Vcount.copy()
    Acell_noisy = mutation_mat.copy()

    missing_rate = args.d
    n_entries = ncells * ncharacters
    nmissing = math.floor(missing_rate * n_entries)
    selected_cell_indices = np.random.randint(ncells, size=nmissing)
    selected_character_indices = np.random.randint(ncharacters, size=nmissing)
    Acell_missing[selected_cell_indices, selected_character_indices] = -1
    Rtotal_missing[selected_cell_indices, selected_character_indices] = 0
    Vcount_missing[selected_cell_indices, selected_character_indices] = 0
    Acell_noisy[selected_cell_indices, selected_character_indices] = -1
    
    # write ground truth files
    prefix = args.o
    with open(f'{prefix}/gt_tree_edgelist.csv', 'w') as out:
        for edge in celltree.edges:
            out.write(f'{edge[0]},{edge[1]}\n')

    with open(f'{prefix}/gt_tree.newick', 'w') as out:
        out.write(tree_to_newick(celltree) + ';')
    
    writeDOT(celltree, f'{prefix}/gt_tree.dot')

    df_B = pd.DataFrame(B, index=list(T.nodes),
                        columns = [f'c{idx}' for idx in range(ncharacters)] + ['cluster_id'], dtype=int)            
    df_Bcell = pd.DataFrame(Bcell, index=[f's{idx}' for idx in range(ncells)],
                            columns = [f'c{idx}' for idx in range(ncharacters)] + ['cluster_id'], dtype=int)            
    df_Acell = pd.DataFrame(Acell, index=[f's{idx}' for idx in range(ncells)],
                            columns = [f'c{idx}' for idx in range(ncharacters)] + ['cluster_id'], dtype=int)    
    df_Acell_noisy = pd.DataFrame(Acell_noisy, index=[f's{idx}' for idx in range(ncells)],
                                  columns = [f'c{idx}' for idx in range(ncharacters)] + ['cluster_id'], dtype=int)
    
    df_Rtotal = pd.DataFrame(Rtotal, index=[f's{idx}' for idx in range(ncells)],
                        columns = [f'c{idx}' for idx in range(ncharacters)], dtype=int)
    df_Vcount = pd.DataFrame(Vcount, index=[f's{idx}' for idx in range(ncells)],
                        columns = [f'c{idx}' for idx in range(ncharacters)], dtype=int)    
    df_Rtotal_missing = pd.DataFrame(Rtotal_missing, index=[f's{idx}' for idx in range(ncells)],
                        columns = [f'c{idx}' for idx in range(ncharacters)], dtype=int)
    df_Vcount_missing = pd.DataFrame(Vcount_missing, index=[f's{idx}' for idx in range(ncells)],
                        columns = [f'c{idx}' for idx in range(ncharacters)], dtype=int)    
    

    df_B.to_csv(f'{prefix}/gt_multi_state_tree_node_character_matrix.csv')
    df_Bcell.to_csv(f'{prefix}/gt_multi_state_character_matrix.csv')
    df_Acell.to_csv(f'{prefix}/gt_character_matrix_without_noise.csv')
    df_Acell_noisy.to_csv(f'{prefix}/gt_character_matrix.csv')
    
    df_Rtotal.to_csv(f'{prefix}/gt_read_count_without_missing.csv')
    df_Vcount.to_csv(f'{prefix}/gt_variant_count_without_missing.csv')
    df_Rtotal_missing.to_csv(f'{prefix}/gt_read_count.csv')
    df_Vcount_missing.to_csv(f'{prefix}/gt_variant_count.csv')
        
    # scarlet data format
    scarlet_readcount_data = []
    scarlet_readcount_data_columns = ['cell_id', 'c']
    for mutation in df_Rtotal_missing.columns:
        scarlet_readcount_data_columns.append(f'{mutation}_v')
        scarlet_readcount_data_columns.append(f'{mutation}_t')

    for cell in df_Rtotal_missing.index:
        cell_readcount_data = [cell, df_Acell.loc[cell]['cluster_id']]
        for mutation in df_Rtotal_missing.columns:
            cell_readcount_data.append(df_Vcount_missing.loc[cell][mutation])
            cell_readcount_data.append(df_Rtotal_missing.loc[cell][mutation])
        scarlet_readcount_data.append(cell_readcount_data)    

    df_scarlet_readcounts = pd.DataFrame(scarlet_readcount_data, columns = scarlet_readcount_data_columns)
    df_scarlet_readcounts = df_scarlet_readcounts.set_index('cell_id')
    df_scarlet_readcounts.to_csv(f'{prefix}/gt_readcounts_scarlet.csv')
    
    with open(f'{prefix}/gt_copynumbertree_scarlet.txt', 'w') as out:    
        for edge in Tc.edges:
            edge_data = [edge[0], edge[1]]
            for mutation in range(nmutations):
                if R[edge[0],mutation] > R[edge[1],mutation]:
                    edge_data.append(f'c{mutation}')
            out.write(','.join(map(str,edge_data)) + '\n')

    # sphyr data format
    with open(f'{prefix}/gt_sphyr.txt', 'w') as out:
        out.write(f'{len(df_Acell_noisy)}\n')
        out.write(f'{len(df_Acell_noisy.columns) - 1}\n')

        for idx, row in df_Acell_noisy.iterrows():
            out.write(f"{' '.join(map(str,row.values[:-1]))}\n")    

    
    # SCITE data format
    Acell_noisy_scite = Acell_noisy[:,:-1].T.copy().astype(int)
    Acell_noisy_scite[Acell_noisy_scite == -1] = 3
    np.savetxt(f"{prefix}/gt_scite.txt", Acell_noisy_scite, delimiter=" ", fmt='%d')

    # SIFIT data format
    np.savetxt(f"{prefix}/gt_sifit.txt", np.hstack((np.arange(Acell_noisy_scite.shape[0])[:,None], Acell_noisy_scite)), delimiter=" ", fmt='%d')    

    # PhiSCS data format
    Acell_noisy_scite_str = Acell_noisy_scite.astype(str)
    Acell_noisy_scite_str = np.char.replace(Acell_noisy_scite_str, '3', '?').T
    snv_mat_phiscs = np.vstack((np.hstack((np.array([['cell_idx/mut_idx']]), np.array(df_Acell_noisy.columns)[np.newaxis,:-1])),
                                np.hstack((np.array(df_Acell_noisy.index)[:,np.newaxis], Acell_noisy_scite_str))))
    np.savetxt(f"{prefix}/gt_phiscs.txt", snv_mat_phiscs, delimiter="\t", fmt='%s')
    
    # Dollo-CDP data format
    Acell_loc = pd.DataFrame(Acell, index=[f's{idx}' for idx in range(ncells)],columns = [f'{names[idx]}' for idx in range(1,ncharacters+1)] + ['cluster_id'], dtype=int)
    write_nexus_file(Acell_loc, f'{prefix}/gt_nexus.nex')
    write_bed_file(Acell_loc, f'{prefix}/gt_locations.bed')


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', type=int, help='number of samples [5]', default = 5)
    parser.add_argument('-m', type=int, help='number of SNV mutations [5]', default = 5)
    parser.add_argument('-p', type=int, help='number of clusters [1]', default = 1)
    parser.add_argument('-k', type=int, help='number of SNV losses per character [0]', default = 0)
    parser.add_argument('-o', type=str, help='output prefix', default='sample')
    parser.add_argument('-s', type=int, help='seed [0]', default = 0)
    parser.add_argument('-d', type=float, help='missing data rate [0.0]', default=0)
    parser.add_argument('--cov', type=int, help='coverage of read count [50]', default = 50)
    parser.add_argument('-a', type=float, help='false positive error rate', default = 0)
    parser.add_argument('-b', type=float, help='false negative error rate', default = 0)
    parser.add_argument('--ado', type=float, help='precision parameter for ado [15]', default = 15)
    parser.add_argument('--maxcn', type=float, help='maximum allowed copy number [8]', default = 8)
    parser.add_argument('--readthreshold', type=int, help='variant read count threshold for generating the mutation matrix [5]', default=5)
    parser.add_argument('--vafthreshold', type=float, help='VAF threshold for generating the mutation matrix [0.1]', default = 0.1)
    parser.add_argument('-l', type=float, help='rate of mutation loss [0.8]', default = 0.8)
    parser.add_argument('-v', action='store_true', default=False)
    args = parser.parse_args(None if sys.argv[1:] else ['-h'])
    
    main(args)
