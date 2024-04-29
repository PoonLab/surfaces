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
        '--prune',  type=str, default=False, 
        help = '<optional> Target length for prunning the trees'
    )

    return parser.parse_args()

def prune_length(phy, target):
    """
    From PoonLab/Bioplus/prunetree.py
    Progressively remove the shortest terminal branches in the tree until
    we reach a target tree length.
    :param phy:  Bio.Phylo object
    :param target:  float, target tree length to prune to
    :return: Bio.Phylo object after prunning 
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

def run_selection_pipeline(alignment, tree_file, cleaned_file, out_name):
    """
    :param alignment: fasta file, nucleotide alignment of a CDS
    :param tree_file: newick file built from the alignment file
    :param cleaned_file: str, name for the cleaned file
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
    seqs = {rec.name: rec.seq for rec in SeqIO.parse(args.ff, "fasta")}
    reader = csv.DictReader(args.clusters)
    clus_seqs = {}  # clustered sequences
            
    # Read clustering results
    for row in reader:
        cluster = row['clusters']
        # Plain name: we need to keep entire header
        # Fastree keeps header before first space
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
        # Label for files related to a cluster
        cluster_label = f"{args.label}_{cluster}"
        
        # Create an independent file with sequences in the cluster
        af_name = f"{cluster_label}.fa"
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

        # Get tree info
        tip_names = set([tip.name for tip in phy.get_terminals()])
        tlen = phy.total_branch_length()
        
        # Prune tree 
        if args.prune:
            target = float(tlen)
            if target <= tlen:
                sys.stderr.write(f"Starting tree length: {tlen}\n")
                pruned = prune_length(phy, target=target)
                af_pruned_name = f"{cluster_label}.pruned.fa"
                
                # write resulting sequences to output file
                for tip in pruned.get_terminals():
                    print(tip.name)
                    record = seqs[tip.name]
                    # _ = SeqIO.write(record, af_pruned_name, "fasta")
            else: 
                sys.stderr.write(f"ERROR: Target length ({target}) \
                                 is smaller than tree length ({tlen})\n")
                pruned = phy
                af_pruned_name = af_name
        
        Phylo.write(pruned, f"{cluster_label}.tree", 'newick')
        
        # Measure selection with FUBAR
        if args.run_sel:
            run_selection_pipeline( af_pruned_name,
                                    f"{cluster_label}.tree",
                                    f"{cluster_label}.cleaned.mafft.fa", 
                                    f"{cluster_label}.FUBAR.json")
        
        sys.exit()
        


