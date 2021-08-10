"""
Read in Newick files and calculate the tree length 
treelength: sum of branch lengths in the tree
"""

import os 
import re
from Bio import Phylo

#dir = "/home/sareh/data/fasttree_output"
dir = "/home/sareh/surfaces/find_cds/working/tree_for_fasttree"
outfile_path = "/home/sareh/data/neighbor_cds/treelength2_cds.csv"
outfile = open(outfile_path,"w+")

for file_name in os.listdir(dir):
    name = file_name.replace("cutter_cds_","")
    
    file_path =os.path.join(dir,file_name)
    
    # read in tree from file
    tr = Phylo.read(file_path, 'newick')
    sys.setrecursionlimit(100000) #idk what this is

    tips = tr.get_terminals()

    total=0
    for tip in tips:
        if tip.branch_length is None:
            continue
        total += tip.branch_length
        print(tip.branch_length)

    outfile.write("{},{},{}\n".format(name,total,len(tips)))
    print("total is {}".format(total))
    
