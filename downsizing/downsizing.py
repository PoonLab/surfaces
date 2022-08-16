##NEED Debugging
"""
input: all nbr genome
output : pruned_nbr_genome & pruned_nbr_accn
1. nbr_dir
2. ref_dir
3. mini_out_dir
4. ft_out_dir
5. pruned_ft_dir
6. pruned_fa_dir
7. pruned_nbr_accn.txt
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

mm2bin = "/home/sareh/2020/scripts/downsizing/minimap2.py"

# INPUT
nbr_dir = "/home/sareh/2020/july/all_nbr_genome"
#nbr_dir = "/home/sareh/2020/sequences/nuc_genome/"
ref_dir = "/home/sareh/2020/july/poly_ref_genome"
#ref_dir = "/home/sareh/2020/sequences/ref_nuc_genome/"
# OUTPUT
base_dir = "/home/sareh/2020/july"
mini_out_dir = "msa_nbr_genome"
#mini_out_dir = "/home/sareh/2020/downsizing/all_genome_msa/"
ft_out_dir = "tr_nbr_genome"
#ft_out_dir = "/home/sareh/2020/downsizing/all_genome_fasttree/"
pruned_ft_dir = "pruned_tr_nbr_genome"
#pruned_ft_dir = "/home/sareh/2020/downsizing/pruned_fasttree"
pruned_fa_dir = "pruned_nbr_genome"
#pruned_fa_dir = "/home/sareh/2020/downsizing/pruned_ncbi_genome"
pruned_accn_dir = "pruned_nbr_accn"

# mpi
from mpi4py import MPI
the_world = MPI.COMM_WORLD
my_number = the_world.Get_rank()
total_number = the_world.Get_size()

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

    tips = tr.get_terminals()

    if target >= len(tips):
        print("Tree already has fewer tips than target.")
        sys.exit()

    while len(tips) > target:
        # find shortest tip
        lengths = [(tip.branch_length, tip.name) for tip in tips]
        lengths.sort()
        min_length, min_tip = lengths[0]
        tr.prune(min_tip)
        tips = tr.get_terminals()

    return(tr)

def accn_from_tre (file):

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
    count = 0
    for ref_file in os.listdir(ref_dir):
        count += 1 # mpi
        if count % total_number != my_number:
                continue
        #file paths
        ref_path = os.path.join(base_dir, ref_dir, ref_file)
        nbr_path = os.path.join(base_dir, nbr_dir, ref_file) # NC_009996.csv
        #accn_fa = "{}.fa".format(ref_file)
        mini_out_path = os.path.join(base_dir, mini_out_dir, ref_file)
        ft_out_path = os.path.join(base_dir, ft_out_dir, ref_file)
        pruned_ft_path = os.path.join(base_dir, pruned_ft_dir, ref_file)
        pruned_fa_path = os.path.join(base_dir, pruned_fa_dir, ref_file)
        pruned_accn_path = os.path.join(base_dir, pruned_accn_dir, ref_file)

        print(mini_out_path)
        nthread = 5

        #1. python3 minimap2.py -o minimap_alignment.fasta --ref ref.hiv.test.fa NC_001802.txt -f -a
        cmd_minimap = ['python3', mm2bin, '-o', str(mini_out_path), '--ref', str(ref_path), str(nbr_path), '-f', '-a', '-t', str(nthread)]
        subprocess.call(cmd_minimap)
        print(mini_out_path)

        #2. fasttree -nt -gtr -gamma minimap_alignment.fasta > fasttree.tre
        cmd_fasttree = ['fasttree', '-nt', '-gtr', '-gamma','-quiet', str(mini_out_path)]
        ft = subprocess.run(cmd_fasttree, stdout=subprocess.PIPE)

        #writing it out
        ft_handle = open(ft_out_path, "w")
        ftout = str(ft.stdout)[2:-3] #remove: b' #remove: '\n
        ft_handle.write(ftout)
        print("fasttree ran!!")
        sleep(2)

        #3. fasttree.tre > pruned_fasttree.tre
        pruned_tr = prune(ft_out_path, 100)
        Phylo.write(pruned_tr, pruned_ft_path, 'newick')
        sleep(2)

        #4. pruned_fasttree.tre > pruned_ncbi_nuc_genome
        pruned_accns = tuple(accn_from_tre(pruned_ft_path))
        sleep(2)

        # write out accns
        pruned_fa = open(pruned_fa_path,'w')
        fasta_file = open(nbr_path, 'r')

        for record in SeqIO.parse(fasta_file, 'fasta'):
            id = []
            x = re.findall('[A-Z]{1,3}[0-9]{5,6}.[0-9]', record.id)

            for i in x:
                if i in accns:
                    SeqIO.write(record,pruned_fa,'fasta')

if __name__ == "__main__":
    main()
    
    
    
    
