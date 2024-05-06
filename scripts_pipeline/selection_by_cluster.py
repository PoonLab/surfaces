import csv
import argparse
import sys
import subprocess
import tempfile
import os

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

def align_amino(records, bin='mafft'):
    """
    Function from Bioplus/codon_align.py
    """
    handle = tempfile.NamedTemporaryFile(delete=False)
    for record in records:
        aaseq = record.seq.translate()
        line = f">{record.name}\n{aaseq}\n"
        handle.write(line.encode('utf-8'))
    handle.close()

    stdout = subprocess.check_output([bin, '--quiet', handle.name])
    stdout = stdout.decode('utf-8')
    os.remove(handle.name)  # delete temp file
    return SeqIO.parse(StringIO(stdout), "fasta")

def apply_aln(records, aln):
    """
    Function from Bioplus/codon_align.py
    Apply the alignment of translated amino acid sequences to the nucleotide
    sequences to yield a codon-aware alignment.
    :param records:  list, SeqRecord objects representing unaligned nucleotide
                     sequences
    :param aln:  SeqIO iterable object from align_amino()
    :yield:  tuple, sequence name and codon-aligned nucleotide sequence
    """
    for i, aaseq in enumerate(aln):
        nucseq = records[i].seq  # original nucleotide sequence
        newseq = ''
        pos = 0
        for aa in aaseq:
            if aa == '-':
                newseq += '---'
                continue
            newseq += nucseq[pos:(pos+3)]
            pos += 3
        yield records[i].name, newseq

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
        # print(f"\n>>>> len: {tlen}, number of tips: {len(tips)}")
        # we have to re-sort every time because removing a branch
        # will lengthen another branch
        tips = sorted(tips, key=lambda x: x.branch_length)
        _ = phy.prune(tips[0])
        tlen -= tips[0].branch_length
        tips = tips[1:]  # update list
    return phy

def run_fasttree(aln):
    # Create a tree from cleaned a nucleotide alignment
    p = subprocess.Popen(['fasttree', '-nt', aln], 
                            stdout=subprocess.PIPE)
    
    tree, err = p.communicate()

    # tree = subprocess.run(['fasttree', '-nt'], stdout=subprocess.PIPE, 
    #                         stderr=subprocess.DEVNULL, stdin=aln)
    # # t = open(tree_file, 'w')
    # t.write(tree.stdout.decode("utf-8"))
    # t.close()
    return(tree)

def align_codon_aware(cluster_seqs):
    """
    Align sequences that belong to the same cluster
    :param seqs: dict, key by sequence label, values are nt sequences
    :return codon_cluster_aln: StringIO object with aligned sequences 
    """
    seqs_string = ""  
    for header, seq in cluster_seqs.items():
        line = f">{header}\n{seq}\n"
        seqs_string += line

    # load records into memory because we need to iterate twice
    seqs_cluster = list(SeqIO.parse(StringIO(seqs_string), "fasta"))
    cluster_aln = align_amino(seqs_cluster)
    
    # Store output as StrinIO object
    codon_cluster_aln = StringIO()
    for label, codseq in apply_aln(seqs_cluster, cluster_aln):
        codon_cluster_aln.write(f">{label}\n{codseq}\n")

    return(codon_cluster_aln)

def run_mafft(alignment):
    # Align nucleotide sequences
    aligned = subprocess.Popen(['mafft', '--quiet', alignment], 
                               stdout=subprocess.PIPE)
    # stdout = stdout.decode('utf-8')
    return(aligned)
    # cmd = f"mafft --quiet {alignment} > {alignment}.mafft"
    # subprocess.call(cmd, shell=True)

def run_selection_pipeline(alignment, label):
    """
    :param alignment: str nucleotide alignment of a CDS
    :param label: str, label for the output files
    """

    # Clean stop codons from alignment
    clean_aln = subprocess.run(['hyphy', 'cln', 'Universal', alignment, 
                                'No/Yes', f"{label}.hyphy.cleaned.fasta"], 
                                stdout=subprocess.PIPE)

    # Re create tree: hyphy modify aln headers, sequences might not match tree anymore
    tree = subprocess.run(['fasttree', '-nt', f"{label}.hyphy.cleaned.fasta"], 
                            stdout=subprocess.PIPE, 
                            stderr=subprocess.DEVNULL)
    
    # Write tree file
    with open(f"{label}.hyphy.tree", 'w') as t_f:
        t_f.write(tree.stdout.decode("utf-8"))
    
    # Run FUBAR
    print(f"--- Running FUBAR for: {label} ---")
    p = subprocess.run(['hyphy', 'fubar', f"{label}.hyphy.cleaned.fasta",
                         f"{label}.hyphy.tree", '--output', f"{label}.FUBAR.json"])
    
    # TO DO: How can we skip generating all this files!!

if __name__=="__main__":
    args = parse_args()
    seqs = {rec.name: rec.seq for rec in SeqIO.parse(args.ff, "fasta")}
    n_prots = args.n_prots
            
    # Read clustering results
    reader = csv.DictReader(args.clusters)
    # Separate sequences from clustering results
    all_clus_seqs = {} 
    for row in reader:
        cluster = row['clusters']
        label = row['name']  
        if n_prots:
            if int(cluster) > n_prots:
                print(f" Skipping {label} from cluster: {cluster}")
                continue

        if cluster not in all_clus_seqs:
            all_clus_seqs[cluster] = {}
        all_clus_seqs[cluster][label] = seqs[label]

    bad_trees = []  # Clusters with tree errors    
    # Align, make tree, prune tree, measure selection
    for cluster in all_clus_seqs:
        
        print(f"\n>>>> Processing cluster: {cluster} <<<<\n")
        clust_seqs = all_clus_seqs[cluster]
        
        # Codon-aware sequence alignment
        codon_cluster_aln = align_codon_aware(clust_seqs)
        
        # Store alignment file so fasttree can process it
        cluster_label = f"{args.label}_{cluster}"
        aln_name = f"{cluster_label}.codon_aln.mafft.fasta"
        with open(aln_name, 'w') as file:
            file.write(codon_cluster_aln.getvalue())

    
        with open(f"{cluster_label}_before_prun.fasta", 'w') as file:
            file.write(codon_cluster_aln.getvalue())
        # Build tree
        tree = run_fasttree(aln_name)
        handle = StringIO(tree.decode())
        phy = Phylo.read(handle, "newick")
        Phylo.write(phy, f"{cluster_label}.before_prun.tree", 'newick')

        # Get tree info
        tip_names = set([tip.name for tip in phy.get_terminals()])
        tlen = phy.total_branch_length()

        # Prune tree 
        if args.prune:
            target = float(args.prune)
            if target < tlen:

                sys.stderr.write(f"\nStarting tree length: {tlen}\n")

                try:
                    # Prune tree
                    pruned_tree = prune_length(phy, target=target)
    
                except Exception as e:
                    print(f"\n-----------------------------------------------")
                    print(f"\tError while prunning tree of cluster {cluster}:")
                    print(f"\t'{e}'")
                    print(f"\tConsider removing sequences from longest branches")
                    print(f"-------------------------------------------------\n")
                    bad_trees.append(cluster)
                    continue
                
                # Pruned tree is now too short
                if pruned_tree.total_branch_length() < (target - 0.05):
                    new_l = pruned_tree.total_branch_length()
                    n_seqs = len(pruned_tree.get_terminals())
                    print(f"\n-----------------------------------------------")
                    print(f"\tError after prunning tree of cluster {cluster}")
                    print(f"\tPruned tree length = {new_l}, number of sequences = {n_seqs}")
                    print(f"\tConsider removing sequences from longest branches")
                    print(f"-------------------------------------------------\n")
                    bad_trees.append(cluster)
                    continue

                # Get sequences after prunning
                pruned_seqs = {}
                for tip in pruned_tree.get_terminals():
                    pruned_seqs[tip.name] = clust_seqs[tip.name]
                
                # Re-align sequences on tips
                pruned_aln = align_codon_aware(pruned_seqs)

                # Overwrite alignment with pruned_tree sequences
                with open(aln_name, 'w') as file:
                    file.write(pruned_aln.getvalue())
                    
            else: 
                print(f"\n>> Target length ({target}) is shorter than tree length ({round(tlen, 3)})\n")       
                pruned_tree = phy
        
        else: 
            pruned_tree = phy

        print(f"\n-------------------------------------------------")
        print(f"\tFinal tree length: {pruned_tree.total_branch_length()}")
        print(f"\tNumber of sequences: {len(pruned_tree.get_terminals())}")
        print(f"\tTree file at: {cluster_label}.tree")
        print(f"-------------------------------------------------\n")
        Phylo.write(pruned_tree, f"{cluster_label}.tree", 'newick')
        
        # Measure selection with FUBAR
        if args.run_sel:
            run_selection_pipeline(aln_name,
                                   cluster_label)
        
        # At the end of for loop, delete alignment from memory
        codon_cluster_aln.close
    
    print(f"\n>>> Unsuccesful analysis for clusters: {bad_trees}\n")