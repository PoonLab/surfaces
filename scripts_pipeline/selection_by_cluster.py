"""
Create multiple fasta, one for each cluster of proteins
"""
import csv
import argparse
import sys
import subprocess
from Bio import SeqIO

def parse_args():
    parser = argparse.ArgumentParser(
        description="Create multiple fasta, one for each cluster of proteins"
    )
    parser.add_argument(
        'clusters', type=argparse.FileType('r'),
        help='Path to file containing cluster output in csv format'
    )
    parser.add_argument(
        'ff', default = None, help = 'Path to fasta file with CDSs'
    )
    parser.add_argument(
        '--label', type=str, default='cluster', help = 'label for CDSs files'
    )
    parser.add_argument(
        '--n_prots', type=int, default=10, 
        help = 'Maximum numer of clusters allowed'
    )
    parser.add_argument(
        '--run_sel', action='store_true',
        help = 'Run selection pipeline for all clustered sequences'
    )

    return parser.parse_args()


def run_selection_pipeline(alignment, cleaned_file, tree_file, out_name):
    """
    :param alignment: file, nucleotide alignment of a CDS
    :param cleaned_file: str, name for the cleaned file
    :param tree_file: str, name of file for CDS tree
    :out_name: str, name and location for FUBAR json file 
    """
    # Align nucleotide sequences
    cmd = f"mafft --quiet {alignment} > {alignment}.mafft"
    subprocess.call(cmd, shell=True)

    # Clean stop codons from alignment
    subprocess.run(['hyphy', 'cln', 'Universal', f"{alignment}.mafft", 'No/Yes', cleaned_file],
    stdout=subprocess.DEVNULL)

    # Create a tree from cleaned CDS alignment
    tree = subprocess.run(['fasttree', '-nt', cleaned_file], stdout=subprocess.PIPE, 
    stderr=subprocess.DEVNULL)
    t = open(tree_file, 'w')
    t.write(tree.stdout.decode("utf-8"))
    t.close()
    
    # Run FUBAR
    print(f">>> Running FUBAR for: {out_name}")
    p = subprocess.run(['hyphy', 'fubar', cleaned_file, tree_file, '--output', out_name],
    stdout=subprocess.DEVNULL)


args = parse_args()
seqs = {rec.description: rec.seq for rec in SeqIO.parse(args.ff, "fasta")}
reader = csv.DictReader(args.clusters)
clus_seqs = {}  # clustered sequences

# Read clustering results
for row in reader:
    cluster = row['clusters']
    label = row['name']
    if cluster not in clus_seqs:
        clus_seqs[cluster] = []
    clus_seqs[cluster].append(label)

# Divide CDSs into different files based on clustering
for cluster in clus_seqs:
    if int(cluster) > args.n_prots:
        continue

    print(f"\nProcesing cluster: {cluster}")
    f_name = f"{args.label}_{cluster}"
    file = open(f"{f_name}.fa", "w")
    for name in clus_seqs[cluster]:
        file.write(f'>{name}\n{seqs[name]}\n')
    file.close()

    # Measure selection with FUBAR
    if args.run_sel:
        run_selection_pipeline( f"{f_name}.fa",
                                f"{f_name}.cleaned.fa",
                                f"{f_name}.cleaned.tree", 
                                f"{f_name}.FUBAR.json")

