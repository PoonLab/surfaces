"""
#! /usr/bin/python
read newick file to get the list of accns after pruning
write out a fasta file with the pruned accn sequences
python3 prune_fasta.py --tree --fasta --outfile
"""

import re
import argparse
from Bio import SeqIO

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--tree', type=str,
                        help='pruned phylogenetic tree in newick format')
    parser.add_argument('--fasta', type=str,
                        help='pre-pruning fasta file')
    parser.add_argument('--outfile', type=str,
                        help='output fasta file')
    return parser.parse_args()

def prune_accn(file):

    """
    takes in path to a newick file
    return list of accns of pruned sequences
    """

    accns = []
    with open(file, 'r') as f:
        all_f = f.read()
        x = re.findall('[A-Z]{1,3}[0-9]{5,6}.[0-9]', all_f)
        accns += x
    return(accns)

def main():
    args = parse_args()
    accns = tuple(prune_accn(args.tree))

    # Opening the files
    outfile = open(args.outfile,'w')
    fasta_file = open(args.fasta, 'r') 
 
    for record in SeqIO.parse(fasta_file, 'fasta'):
        id = []
        x = re.findall('[A-Z]{1,3}[0-9]{5,6}.[0-9]', record.id)

        for i in x:
            if i in accns:
                SeqIO.write(record,outfile,'fasta')
   
    #closing the file  
    outfile.close()
    fasta_file.close()
if __name__ == "__main__":
    main()
