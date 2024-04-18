import csv
import argparse
import sys
import subprocess
from Bio import SeqIO, Phylo
from io import StringIO

description = """\
Create separate fasta file based on cluster classification.
Align, create phylogeny, and measure selection with FUBAR. 
Prune tree to a given length when specified. 
"""

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
        '--n_prots', type=int, default=False, 
        help = 'Maximum number of clusters allowed'
    )
    parser.add_argument(
        '--run_sel', action='store_true',
        help = 'Run selection pipeline for all clustered sequences'
    )
    parser.add_argument(
        '--pn',  type=str, default=False, 
        help = 'Prune tree'
    )

    return parser.parse_args()

def run_prunetree(phy, aln):
    """
    Function from PoonLab/Bioplus/prunetree.py
    Progressively remove the shortest terminal branches in the tree until
    we reach a target tree length.
    :param phy:  Bio.Phylo object
    :param aln:  Sequence alignment
    """
    pruned = subprocess.run['python3', 'prunetree.py', phy, '--mode treelen', '-t 0.5', f'--seq {aln}', f'--outfile {aln}.pruned']
    return pruned

def prune_length(phy, target):
    """
    From PoonLab/Bioplus/prunetree.py
    Progressively remove the shortest terminal branches in the tree until
    we reach a target tree length.
    :param phy:  Bio.Phylo object
    :param target:  float, target tree length to prune to
    :return:  Bio.Phylo object
    """
    tlen = phy.total_branch_length()
    if target >= tlen:
        sys.stderr.write(f"prune_length: requirement already met "
                         f"({tlen}<={target}")
    tips = phy.get_terminals()
    while tlen > target:
        # we have to re-sort every time because removing a branch
        # will lengthen another branch
        tips = sorted(tips, key=lambda x: x.branch_length)
        _ = phy.prune(tips[0])
        tlen -= tips[0].branch_length
        tips = tips[1:]  # update list
    return phy


def run_fasttree(aln):
    # Create a tree from cleaned a nucleotide alignment
    tree = subprocess.run(['fasttree', '-nt'], stdout=subprocess.PIPE, 
                            stderr=subprocess.DEVNULL, stdin=aln.stdout)
    # t = open(tree_file, 'w')
    # t.write(tree.stdout.decode("utf-8"))
    # t.close()
    return(tree)
    
def run_mafft(alignment):
    # Align nucleotide sequences
    stdout = subprocess.Popen(['mafft', '--quiet', alignment], stdout=subprocess.PIPE)
    # stdout = stdout.decode('utf-8')
    return(stdout)
    # cmd = f"mafft --quiet {alignment} > {alignment}.mafft"
    # subprocess.call(cmd, shell=True)

def run_selection_pipeline(alignment, cleaned_file, tree_file, out_name):
    """
    :param alignment: file, nucleotide alignment of a CDS
    :param cleaned_file: str, name for the cleaned file
    :param tree_file: str, name of file for CDS tree
    :out_name: str, name and location for FUBAR json file 
    """

    # Clean stop codons from alignment
    subprocess.run(['hyphy', 'cln', 'Universal', f"{alignment}.mafft", 'No/Yes', cleaned_file],
    stdout=subprocess.DEVNULL)

    # Run FUBAR
    print(f">>> Running FUBAR for: {out_name}")
    p = subprocess.run(['hyphy', 'fubar', cleaned_file, tree_file, '--output', out_name],
    stdout=subprocess.DEVNULL)

if __name__=="__main__":
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
    
    n_prots = args.n_prots if args.n_prots else len(clus_seqs)

    # Divide CDSs into different files based on clustering
    for cluster in clus_seqs:
        if int(cluster) > n_prots:
            continue

        print(f"\nProcesing cluster: {cluster}")
        af_name = f"{args.label}_{cluster}.fa"
        file = open(af_name, "w")
        for name in clus_seqs[cluster]:
            file.write(f'>{name}\n{seqs[name]}\n')
        file.close()

        # Align sequences
        align = run_mafft(af_name)

        # Build tree
        tree = run_fasttree(align)

        handle = StringIO(tree.stdout.decode("utf-8"))
        phy = Phylo.read(handle, "newick")
        tlen = phy.total_branch_length()
        print(f"Tree length: {phy.total_branch_length()}")
        
        if args.pn:
            target = float(tlen)
            if target <= tlen:
                sys.stderr.write(f"Starting tree length: {tlen}\n")
                sys.exit()
            pruned = prune_length(phy, target=target)
            
        
        
        sys.exit()
        

        # # Measure selection with FUBAR
        # if args.run_sel:
        #     run_selection_pipeline( f"{f_name}.fa",
        #                             f"{f_name}.cleaned.fa",
        #                             f"{f_name}.cleaned.tree", 
        #                             f"{f_name}.FUBAR.json")

