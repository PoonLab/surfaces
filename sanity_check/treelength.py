
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

file_path = sys.argv[1] #/home/sareh/re-do/all_neighbors/fasttree_genome

# read in tree from file
phy = Phylo.read(file_path, 'newick')

total = phy.count_terminals() # Total number of nodes 
tree_length = phy.total_branch_length() # Total branch length

result <- paste(name,total,tree_length)

# Writing it out to a file 
output.file <- file("./gene_treelength.csv", "a")
cat(result, file=output.file, append=TRUE, sep = "\n")
close(output.file)

print("total is {}".format(total))
