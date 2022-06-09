
#!/usr/bin/python3
"""
Read in Newick files and calculate the tree length 
treelength: sum of branch lengths in the tree
biopython total_branch_length()
"""

import sys
import os 
import re
from Bio import Phylo

dir = sys.argv[1] #/home/sareh/re-do/all_neighbors/fasttree_genome
outfile_path = sys.argv[2] #/home/sareh/re-do/treelength/treelength.csv

outfile = open(outfile_path,"w+")

for file_name in os.listdir(dir):
    
    file_path =os.path.join(dir,file_name)
    
    # read in tree from file
    phy = Phylo.read(file_path, 'newick')

    total = phy.count_terminals()   

    tree_length = phy.total_branch_length()

    outfile.write("{},{},{}\n".format(name,total,tree_length))
    print("total is {}".format(total))
