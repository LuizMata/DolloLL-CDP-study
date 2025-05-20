#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 5 2022

@author: Palash Sashittal
"""

import pandas as pd
import sys
import argparse
import itertools
import math
import numpy as np
import dendropy
import os

import networkx as nx
import ete3
from solveConstrainedDollo import solveConstrainedDollo
from dendropy.calculate.treecompare import false_positives_and_negatives

sys.setrecursionlimit(1500)

def compare_trees(tr1, tr2):
    """
    Compares two trees
    Parameters
    ----------
    tr1 : dendropy tree object
            First tree (typically the model tree)
    tr2 : dendropy tree object
            Second tree (typically the estimated tree)
    Returns
    -------
    nl : int
         Size of the shared leaf set, i.e., the number of leaves in both trees
    i1 : int
          Number of internal edges in tree 1 after restriction to shared leaves
    i2 : int
          Number of internal edges in tree 2 after restriction to shared leaves
    fn : int
         Number of edges in tree 1 that are not in tree 2
    fp : int
         Number of edges in tree 2 that are not in tree 1
    rf : float
         Normalized Robinson-Foulds (RF) distance between tree 1 and 2
    Example
    -------
    If tree 1 corresponds to "(((A,B,C),D),E);"
    and tree 2 corresponds to "((((A,B),C),D),E);",
    then the output is "5 1 2 0 1 0.25". In this example,
      + tree 1 and tree 2 share five leaves (A, B, C, D, E).
      + tree 1 has one internal edge "A,B,C|D,E"
      + tree 2 has two internal edges "A,B|C,D,E" and "A,B,C|D,E"
      + no edges in the tree 1 are missing from tree 2
      + one edge in the tree 2 is missing from the tree 1
      + normalized RF distance is (FP+FN)/(2*NL-6) = (1+0)/(2*5-6) = 0.25
    """

    # Unroot the two trees!
    tr1.is_rooted = False
    tr1.collapse_basal_bifurcation(set_as_unrooted_tree=True)

    tr2.is_rooted = False
    tr2.collapse_basal_bifurcation(set_as_unrooted_tree=True)

    # Restrict the two trees to the same leaf set if necessary!
    lb1 = set([l.taxon.label for l in tr1.leaf_nodes()])
    lb2 = set([l.taxon.label for l in tr2.leaf_nodes()])

    com = lb1.intersection(lb2)
    if com != lb1 or com != lb2:
        com = list(com)
        tns = dendropy.TaxonNamespace(com)

        tr1.retain_taxa_with_labels(com)
        tr1.migrate_taxon_namespace(tns)

        tr2.retain_taxa_with_labels(com)
        tr2.migrate_taxon_namespace(tns)
    com = list(com)

    # Compare trees!
    tr1.update_bipartitions()
    tr2.update_bipartitions()

    nl = len(com)
    i1 = len(tr1.internal_edges(exclude_seed_edge=True))
    i2 = len(tr2.internal_edges(exclude_seed_edge=True))

    [fp, fn] = false_positives_and_negatives(tr1, tr2)
    rf = (fn + fp) / (2.0 * nl - 6.0)

    return(nl, i1, i2, fn, fp, rf)

def treescore(tree1string, tree2string):
    taxa = dendropy.TaxonNamespace()

    tree1 = dendropy.Tree.get(path=tree1string,
                              schema='newick',
                              rooting='force-unrooted',
                              taxon_namespace=taxa)

    tree2 = dendropy.Tree.get(data=tree2string,
                              schema='newick',
                              rooting='force-unrooted',
                              taxon_namespace=taxa)

    return compare_trees(tree1, tree2)



def remove_species(newick: str, target: str) -> str:
    def parse_subtree(s):
        stack = []
        current = ''
        for char in s:
            if char == '(':
                stack.append([])
            elif char == ',':
                if current.strip():
                    stack[-1].append(current.strip())
                current = ''
            elif char == ')':
                if current.strip():
                    stack[-1].append(current.strip())
                current = ''
                subtree = stack.pop()
                label = ''
                if stack:
                    stack[-1].append(subtree)
                else:
                    return subtree
            else:
                current += char
        return s.strip(';')  # fallback

    def remove_node(tree):
        if isinstance(tree, str):
            return None if tree == target else tree
        else:
            pruned = [remove_node(child) for child in tree]
            pruned = [child for child in pruned if child is not None]
            if not pruned:
                return None
            elif len(pruned) == 1:
                return pruned[0]
            else:
                return pruned

    def to_newick(tree) -> str:
        if isinstance(tree, str):
            return tree
        else:
            return '(' + ','.join(to_newick(child) for child in tree) + ')'

    parsed = parse_subtree(newick)
    pruned = remove_node(parsed)
    if pruned is None:
        return ';'  # No tree left
    return to_newick(pruned) + ';'

def get_dollocdp_tree_dataframe(dollocdp_newick):

    if not os.path.exists(dollocdp_newick):
        print(f"This job did not finish: {dollocdp_newick}")
        sys.exit(1)

    dollocdp_newick_string = ''
    with open(dollocdp_newick, 'r') as inp:
        for line in inp:
            dollocdp_newick_string += line.rstrip('\n')    

    tree = remove_species(dollocdp_newick_string, "O")

    return tree

def get_scarlet_tree_dataframe(scarlet_fname, scarlet_edgelist_fname):
    
    df_corrected = pd.read_csv(scarlet_fname, index_col = 0)
    
    T_scarlet_nx = nx.DiGraph()
    with open(scarlet_edgelist_fname, 'r') as inp:
        for line in inp:
            nodes = line.rstrip('\n').split(',')
            parent_node = nodes[0].split(' ')[0]
            child_node = nodes[1].split(' ')[0]

            parent_type = parent_node.split(':')[0]
            parent_name = parent_node.split(':')[1]
            child_type = child_node.split(':')[0]
            child_name = child_node.split(':')[1]

            if parent_type == 'ROOT':
                parent_node_name = f'r{parent_name}'
            elif parent_type == 'MUT':
                parent_node_name = f'{parent_name}_1'
            else:
                parent_node_name = parent_name

            if child_type == 'ROOT':
                child_node_name = f'r{child_name}'
            elif child_type == 'MUT':
                child_node_name = f'{child_name}_1'
            else:
                child_node_name = child_name

            T_scarlet_nx.add_edge(parent_node_name, child_node_name)
    
    Tsol = ete3.Tree(tree_to_newick(T_scarlet_nx) + ';')
    
    return Tsol, df_corrected, T_scarlet_nx

def get_condor_tree_dataframe(condor_fname, condor_newick):
    

    if not os.path.exists(condor_newick):
        print(f"This job did not finish: {condor_newick}")
        sys.exit(1)

    condor_newick_string = ""

    with open(condor_newick,'r') as f:
        condor_newick_string += f.read().replace('\n', '')

    
    df_multi = pd.read_csv(condor_fname, index_col = 0)

    df_corrected = df_multi.copy()
    df_corrected[df_corrected > 1] = 0

    Tsol = condor_newick_string
    
    df_binary = solveConstrainedDollo.expand_multi_state_to_binary(df_multi)
    _, nxtree_condor = solveConstrainedDollo.generate_perfect_phylogeny(df_binary)    
    
    return Tsol, df_corrected, nxtree_condor

def read_sphyr(fname):
    
    with open(fname, 'r') as inp:
        idx = 0
        data = []
        for line in inp:
            if idx == 0:
                n = int(line.split(' ')[0])
            elif idx == 1:
                m = int(line.split(' ')[0])
            else:
                data.append(list(map(int, line.split(' '))))
            idx += 1
    
    return pd.DataFrame(data, columns = [f'c{idx}' for idx in range(m)], index = [f's{idx}' for idx in range(n)])

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

def get_sphyr_tree_dataframe(sphyr_fname):

    df_multi = read_sphyr(sphyr_fname)
    df_binary = solveConstrainedDollo.expand_multi_state_to_binary(df_multi)
    _, nxtree_sphyr = solveConstrainedDollo.generate_perfect_phylogeny(df_binary)

    df_corrected = df_multi.copy()
    df_corrected[df_corrected > 1] = 0

    T_sphyr = ete3.Tree(tree_to_newick(nxtree_sphyr)+';', format=1)
    
    return T_sphyr, df_corrected, nxtree_sphyr

def get_sifit_tree_dataframe(sifit_fname, sifit_newick):
    df_corrected = pd.read_csv(sifit_fname, sep='\t', header=None, index_col = 0)
    n = len(df_corrected)
    m = len(df_corrected.columns)
    
    df_corrected.rename(columns={idx: f'c{idx-1}' for idx in range(1,m+1)}, index={f'sc{idx}': f's{idx-1}' for idx in range(1, n+1)}, inplace=True)
    df_corrected = df_corrected.loc[[f's{idx}' for idx in range(n)]]
    df_corrected = df_corrected.sort_index()
    
    sifit_newick_string = ''
    with open(sifit_newick, 'r') as inp:
        for line in inp:
            sifit_newick_string += line.rstrip('\n')    

    mod_sifit_newick_list = []
    for idx, string in enumerate(sifit_newick_string.split('sc')):
        if idx == 0:
            mod_sifit_newick_list.append(string)
        else:
            splited_string = string.split(':')
            splited_string[0] = str(int(splited_string[0]) - 1)
            mod_string = ':'.join(splited_string)
            mod_sifit_newick_list.append('s' + mod_string)    

    T_sifit = ete3.Tree(''.join(mod_sifit_newick_list))
    
    return T_sifit, df_corrected

def get_scite_tree_dataframe(n, m, gv_fname):

    # scite_nx_tree = nx.DiGraph()
    # with open(gv_fname, 'r') as inp:
    #     for line in inp:
    #         if not line.startswith('digraph') and not line.startswith('node') and not line.startswith('}'):
    #             data = line.rstrip(';\n').split(' -> ')
    #             scite_nx_tree.add_edge(data[0], data[1])
    # T_scite = ete3.Tree(tree_to_newick(scite_nx_tree) + ';')
    
    gv_file = open(gv_fname, 'r')
    gv_file.readline()
    gv_file.readline()    

    pi = [-2 for i in range(m+1)]
    for _ in range(m):
        line = gv_file.readline()
        data = line.rstrip(";\n").split()
        source = int(data[0]) - 1
        target = int(data[2]) - 1
        assert 0 <= target < m
        parent = -2
        if source == m:
            parent = -1
        else:
            parent = source
        pi[target] = parent    

    samples = [-2 for i in range(n)]
    gv_file.readline()

    for _ in range(n):
        line = gv_file.readline()
        data = line.rstrip(";\n").split()
        source = int(data[0]) - 1
        target = int(data[2][1:])
        parent = -2
        if source == n:
            parent = -1
        else:
            parent = source
        samples[target] = source        
        
    scite_corrected = []
    for p in range(len(samples)):
        states = [ 0 for c in range(len(pi)) ]
        parent = samples[p]
        while parent != -1:
            #print p, parent
            states[parent] = 1
            parent = pi[parent]
        scite_corrected.append(states[:len(pi) - 1])
        
    df_corrected = pd.DataFrame(scite_corrected, columns = [f'c{idx}' for idx in range(m)], index = [f's{idx}' for idx in range(n)])
    
    df_corrected_modified = df_corrected.copy()
    df_corrected_modified.columns = [f'{x}_1' for x in df_corrected.columns]
    _, nxtree_scite = solveConstrainedDollo.generate_perfect_phylogeny(df_corrected_modified)
    T_scite = ete3.Tree(tree_to_newick(nxtree_scite) + ';')
    
    return T_scite, df_corrected, nxtree_scite

def get_descendant_mutations(T, node):
    if node not in T.nodes:
        return []
    if node.startswith('c'):
        descendants = [node]
    else:
        return []
    for child in T[node]:
        descendants += get_descendant_mutations(T, child)
    return descendants

def get_clustered_mutations(T, node):
    if node not in T.nodes:
        return []
    if node.startswith('c'):
        clustered = [node]
        if len(list(T[node])) == 1:
            child = list(T[node])[0]
            if child.startswith('c'):
                clustered += get_clustered_mutations(T, child)
    else:
        return []
    return clustered

def main(args):

    n = args.n
    m = args.m
    p = args.p
    k = args.k
    d = args.d
    if d == 0:
        d = int(d)
    s = args.s
    
    output_fname = args.o
    
    method = args.method
    simulation_dir = args.dir
    
    if method == 'dollocdp':
        
        #remember to provide this output from DolloCDP
        dollocdp_newick = f'{simulation_dir}/n{n}_m{m}_p{p}_k{k}_d{d}/r{s}/{method}.tre.contracted'

        Tsol = get_dollocdp_tree_dataframe(dollocdp_newick)
        
    elif method == 'dollocdpLL':
        dollocdpll_newick = f'{simulation_dir}/n{n}_m{m}_p{p}_k{k}_d{d}/r{s}/{method}.tre'

        Tsol = get_dollocdp_tree_dataframe(dollocdpll_newick)
        
    elif method == 'condor':
        
        condor_fname = f'{simulation_dir}/n{n}_m{m}_p{p}_k{k}_d{d}/r{s}/{method}_B.csv'
        condor_newick = f'{simulation_dir}/n{n}_m{m}_p{p}_k{k}_d{d}/r{s}/{method}_tree.newick.contracted'
        Tsol, df_corrected, nxtree_sol = get_condor_tree_dataframe(condor_fname, condor_newick)
        
    elif method == 'condorreads':

        condor_fname = f'{simulation_dir}/n{n}_m{m}_p{p}_k{k}_d{d}/r{s}/{method}_B.csv'
        condor_newick = f'{simulation_dir}/n{n}_m{m}_p{p}_k{k}_d{d}/r{s}/{method}_tree.newick'
        Tsol, df_corrected, nxtree_sol = get_condor_tree_dataframe(condor_fname, condor_newick)
        
    elif method == 'sphyr':
        
        sphyr_fname = f'{simulation_dir}/n{n}_m{m}_p{p}_k{k}_d{d}/r{s}/{method}.out'
        Tsol, df_corrected, nxtree_sol = get_sphyr_tree_dataframe(sphyr_fname)

    elif method == 'scarlet':
        
        scarlet_fname = f'{simulation_dir}/n{n}_m{m}_p{p}_k{k}_d{d}/r{s}/{method}.B'
        scarlet_edgelist_fname = f'{simulation_dir}/n{n}_m{m}_p{p}_k{k}_d{d}/r{s}/{method}.edgelist'
        Tsol, df_corrected, nxtree_sol = get_scarlet_tree_dataframe(scarlet_fname, scarlet_edgelist_fname)
    
    # evaluate
    gt_fname = f'{simulation_dir}/n{n}_m{m}_p{p}_k{k}_d{d}/r{s}/gt_character_matrix_without_noise.csv'
    gt_newick = f'{simulation_dir}/n{n}_m{m}_p{p}_k{k}_d{d}/r{s}/gt_tree.newick.contracted'
    gt_multi = f'{simulation_dir}/n{n}_m{m}_p{p}_k{k}_d{d}/r{s}/gt_multi_state_character_matrix.csv'
    df_character_matrix = pd.read_csv(gt_fname, index_col = 0)
    df_character_matrix = df_character_matrix[df_character_matrix.columns[:-1]]
    df_Bcell = pd.read_csv(gt_multi, index_col = 0)
    _, nxtree_gt = solveConstrainedDollo.generate_perfect_phylogeny(solveConstrainedDollo.expand_multi_state_to_binary(df_Bcell[df_Bcell.columns[:-1]]))
    
    #metric
    
    [nl, i1, i2, fn, fp, rf] = treescore(gt_newick, Tsol)

    # mutation error
    nerror = 0
    
    # recall for ancestry and incomparibility
    if (method != 'dollocdp' and method != 'dollocdpLL'):
        sol_descendant_dictionary = {}
        for mutation_idx in range(m):
            sol_descendant_dictionary[f'c{mutation_idx}_1'] = get_descendant_mutations(nxtree_sol, f'c{mutation_idx}_1')

        sol_clustered_dictionary = {}
        for mutation_idx in range(m):
            sol_clustered_dictionary[f'c{mutation_idx}_1'] = get_clustered_mutations(nxtree_sol, f'c{mutation_idx}_1')            
            
        gt_descendant_dictionary = {}
        for mutation_idx in range(m):
            gt_descendant_dictionary[f'c{mutation_idx}_1'] = get_descendant_mutations(nxtree_gt, f'c{mutation_idx}_1')    

        gt_clustered_dictionary = {}
        for mutation_idx in range(m):
            gt_clustered_dictionary[f'c{mutation_idx}_1'] = get_clustered_mutations(nxtree_gt, f'c{mutation_idx}_1')
            
        confusion_mat = np.zeros((4,4))
        for mutation_idx1, mutation_idx2 in itertools.combinations(range(m), 2):
            mutation1 = f'c{mutation_idx1}_1'
            mutation2 = f'c{mutation_idx2}_1'

            if mutation2 in sol_descendant_dictionary[mutation1] and mutation1 in sol_descendant_dictionary[mutation2]:
                print('problem sol')
                print(mutation1, mutation2)
                break    

            x_idx = 3
            if mutation1 in sol_clustered_dictionary[mutation2] or mutation2 in sol_clustered_dictionary[mutation1]:
                x_idx = 2
            else:
                if mutation2 in sol_descendant_dictionary[mutation1]:
                    x_idx = 0
                if mutation1 in sol_descendant_dictionary[mutation2]:
                    x_idx = 1

            if mutation2 in gt_descendant_dictionary[mutation1] and mutation1 in gt_descendant_dictionary[mutation2]:
                print('problem gt')
                print(mutation1, mutation2)
                break

            y_idx = 3
            if mutation1 in gt_clustered_dictionary[mutation2] or mutation2 in gt_clustered_dictionary[mutation1]:
                y_idx = 2
            else:
                if mutation2 in gt_descendant_dictionary[mutation1]:
                    y_idx = 0
                if mutation1 in gt_descendant_dictionary[mutation2]:
                    y_idx = 1

            confusion_mat[x_idx, y_idx] += 1

        ancestry_recall = (confusion_mat[0,0] + confusion_mat[1,1]) /  np.sum(confusion_mat[:,:2])
        # incomparability_recall = confusion_mat[2,2] / np.sum(confusion_mat[:,2])
        incomparability_recall = np.nan
        accuracy = np.trace(confusion_mat)/np.sum(confusion_mat)                
    else:
        ancestry_recall = np.nan
        incomparability_recall = np.nan
        accuracy = np.nan
    
    tp = 0
    fnr = 0
    fpr = 0
    precision = 0
    recall = 0

    if i1 == 0 and i2 == 0:
        tp = nl
    elif i1 == 0 and i2 != 0:
        tp = nl
        fpr = fp/i2
    elif i1 != 0 and i2 == 0:
        tp = nl
        fnr = fn/i1
        fpr = 0
    elif i1 != 0 and i2 != 0:
        tp = (i1-fn) + nl
        fnr = fn/i1
        fpr = fp/i2

    precision = tp/(tp+fp)
    recall = tp/(tp+fn)


    df_result = pd.DataFrame([[n, m, p, k, d, s, method, fp, fn, fpr, fnr, precision, recall]],
                             columns = ['ncells', 'ncharacters', 'nclusters', 'k', 'dropout', 'seed', 'method', 'FP', 'FN', 'FPR', 'FNR', 'Precision', 'Recall'])
    
    df_result.to_csv(f'{output_fname}')
    
def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', type=int, help='number of samples [5]', default = 5)
    parser.add_argument('-m', type=int, help='number of SNV mutations [5]', default = 5)
    parser.add_argument('-p', type=int, help='number of clusters [1]', default = 1)
    parser.add_argument('-k', type=int, help='number of SNV losses per character [0]', default = 0)
    parser.add_argument('-o', type=str, help='output filename', default='sample.csv')
    parser.add_argument('-s', type=int, help='seed [0]', default = 0)
    parser.add_argument('-d', type=float, help='missing data rate [0.0]', default=0)
    parser.add_argument('--dir', type=str, help='simulation directory', required=True)
    parser.add_argument('--method', type=str, help='method for comparison with ground truth', required=True)
    parser.add_argument('-v', action='store_true', default=False)
    args = parser.parse_args(None if sys.argv[1:] else ['-h'])

    main(args)
