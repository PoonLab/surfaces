#! /usr/bin/python

import sys
from Bio import Phylo

phy = Phylo.read(sys.argv[1], 'newick')
thrshld = float(sys.argv[2]) # normalized threshold 
outpath = sys.argv[3]

tips = phy.get_terminals()

tbl = phy.total_branch_length()
genomes = phy.count_terminals()
norm_tbl = tbl/genomes

if thrshld > norm_tbl:
    print("{} Treelength is less than threshold".format(sys.argv[1]))
    sys.exit()

for branch_length, tip_name in lengths:
    print(tip_name)
    lengths = [(tip.branch_length, tip.name) for tip in tips] 
    # [(branchlength, Accession),...]
    lengths.sort() # ascending order 
    phy.prune(tip_name)
    tbl -= branch_length # update (tbl) tree length (total branch length)
    genomes -= 1 # update number of genome
    norm_tbl = tbl/genomes #normalized threshold 
    if norm_tbl < thrshld:
        break

print("{},{},{},this is printing\n".format(sys.argv[1],genomes,tbl))
Phylo.write(phy,outpath, 'newick')
