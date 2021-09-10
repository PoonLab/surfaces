#! /usr/bin/python

import sys
from Bio import Phylo

phy = Phylo.read(sys.argv[1], 'newick')
thrshld = float(sys.argv[2])
outpath = sys.argv[3]

tips = phy.get_terminals()

tree_length = phy.total_branch_length()

if thrshld > tree_length:
    print("{} Treelength is less than threshold".format(sys.argv[1]))
    sys.exit()

lengths = [(tip.branch_length, tip.name) for tip in phy.get_terminals()]
lengths.sort()  # ascending order
for _, tip_name in lengths:
    phy.prune(tip_name)
    if phy.total_branch_length() < thrshld:
        break

print("{},{},{}\n".format(sys.argv[1],phy.count_terminals(),tree_length))
Phylo.write(phy,outpath, 'newick')
