"""
takes in ft output prune it and return both accns and pruned fasta file
INPUT: nbr_dir | fasta files of full nbr genomes
INPUT: ft_out_dir | Newick fasttree output of msa of nbr genome 
OUTPUT: pruned_ft_dir = "pruned_tr_nbr_genome" | Newick
OUTPUT: pruned_accn_dir = "pruned_nbr_accn" | one accn on each line
OUTPUT: pruned_fa_dir = "pruned_nbr_genome" | fasta 
"""

#!/usr/bin python3

import sys
from Bio import Phylo
import subprocess
import os
import re
import argparse
from Bio import SeqIO
import time

# INPUT
base_dir = "/home/sareh/all/checking"
nbr_dir = "/home/sareh/all/nbr_genome"
ft_out_dir = "ft"

# OUTPUT
pruned_ft_dir = "pruned_tr_nbr_genome"
pruned_accn_dir = "pruned_nbr_accn"
pruned_fa_dir = "pruned_nbr_genome"

def convert_fasta (handle):
    """
    takes in handel to an open fasta file
    """
    result = []
    sequence = ''
    # handle = open(file,'r')
    for line in handle:
        if line.startswith('$'): # skip header line
            continue
        elif line.startswith('>') or line.startswith('#'):
            if len(sequence) > 0:
                result.append([h,sequence])
                sequence = ''   # reset
            h = line.strip('>#\n')
        else:
            sequence += line.strip('\n').upper()

    result.append([h,sequence]) # handle last entry
    return result

def prune (tr_path, target):
    tr = Phylo.read(tr_path, 'newick')
    #<class 'Bio.Phylo.Newick.Tree'>
    #{'root': Clade(confidence=1.0), 'rooted': False, 'id': None, 'name': None, 'weight': 1.0}
    #['__class__', '__delattr__', '__dict__', '__dir__', '__doc__', '__eq__', '__format__', '__ge__', '__getattribute__', '__gt__', '__hash__', '__init__', '__init_subclass__', '__le__', '__lt__', '__module__', '__ne__', '__new__', '__reduce__', '__reduce_ex__', '__repr__', '__setattr__', '__sizeof__', '__str__', '__subclasshook__', '__weakref__', '_filter_search', 'as_phyloxml', 'clade', 'collapse', 'collapse_all', 'common_ancestor', 'count_terminals', 'depths', 'distance', 'find_any', 'find_clades', 'find_elements', 'format', 'from_clade', 'get_nonterminals', 'get_path', 'get_terminals', 'id', 'is_bifurcating', 'is_monophyletic', 'is_parent_of', 'is_preterminal', 'is_terminal', 'ladderize', 'name', 'prune', 'randomized', 'root', 'root_at_midpoint', 'root_with_outgroup', 'rooted', 'split', 'total_branch_length', 'trace', 'weight']

    tips = tr.get_terminals()

    if target >= len(tips):
        print("Tree already has fewer tips than target.")
        # break

    while len(tips) > target:
        # find shortest tip
        lengths = [(tip.branch_length, tip.name) for tip in tips]
        lengths.sort()
        min_length, min_tip = lengths[0]
        tr.prune(min_tip)
        tips = tr.get_terminals()
    return(tr)

def main():
    count = 0
    ft_dir = os.path.join(base_dir, ft_out_dir)
    for file in os.listdir(ft_dir):
        print(file)
        #file paths
        nbr_path = os.path.join(base_dir, nbr_dir, file) # NC_009996.csv
        #accn_fa = "{}.fa".format(ref_file)
        ft_out_path = os.path.join(base_dir, ft_out_dir, file)
        pruned_ft_path = os.path.join(base_dir, pruned_ft_dir, file)
        pruned_fa_path = os.path.join(base_dir, pruned_fa_dir, file)
        pruned_accn_path = os.path.join(base_dir, pruned_accn_dir, file)

        nthread = 5

        #prune
        pruned_tr = prune(ft_out_path, 100) #Clade(branch_length=0.000721829, name='KT121671')

        #3. fasttree.tre > pruned_fasttree.tre
        Phylo.write(pruned_tr, pruned_ft_path, 'newick')
        time.sleep(1)

        #4. pruned_fasttree.tre > pruned_ncbi_nuc_genome
        tips = pruned_tr.get_terminals()
        pruned_accns = []
        for tip in tips:
            pruned_accns.append(tip.name)
        #writing it out
        pruned_accns_handle = open(pruned_accn_path, "a")
        for a in pruned_accns:
            pruned_accns_handle.write(a + "\n")

        # nbr fasta
        nbr_handle = open(nbr_path, 'r')
        pruned_fa = open(pruned_fa_path,'w')
        print(nbr_path)

        for record in SeqIO.parse(nbr_handle, 'fasta'):
            # print(record.id) # MK440626,Sindbis
            # x = re.findall('[A-Z]{1,3}[0-9]{5,6}.[0-9]', record.id)
            id = record.id.split(",")[0]
            print(id)
            if id in pruned_accns:
                SeqIO.write(record,pruned_fa,'fasta')

if __name__ == "__main__":
    main()
