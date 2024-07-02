import csv
import argparse
import sys
import subprocess
import tempfile
import os

from Bio import SeqIO, Phylo
from io import StringIO
from glob import glob

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
        'cds_file', nargs='*',
        help = 'Path to fasta file with CDSs'
    )
    parser.add_argument(
        '--clusters_info', '-ci', type=argparse.FileType('r'), default=False,
        help='Path to file containing cluster output in csv format'
    )
    parser.add_argument(
        '--label', '-l' , type=str, default='cluster',
        help = 'label for CDSs files'
    )
    parser.add_argument(
        '--n_prots', type=int, default=False, 
        help = 'Maximum number of clusters allowed'
    )
    parser.add_argument(
        '-s', '--run_sel', action='store_true',
        help = 'Run selection pipeline for all clustered sequences'
    )
    parser.add_argument(
        '-p', '--prune', type=str, default=False,
        help = '<optional> Target length for prunning the trees'
    )
    parser.add_argument(
        '-aln', '--save_before', action='store_true',
        help = 'Save alignment before prunning the tree'
    )
    parser.add_argument(
        '-sp', '--save_prune_res', action='store_true',
        help = 'Store cluster size and tree length before and after prunning in table'
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

def prune_length(phy, target, range):
    """
    From PoonLab/Bioplus/prunetree.py
    Progressively remove the shortest terminal branches in the tree until
    we reach a target tree length.
    :param phy:  Bio.Phylo object
    :param target:  float, target tree length to prune to
    :param range: float, in case users want more flexibility 
                        so the pruned tree could be a bit longer than target
    :return: Bio.Phylo object after pruning 
    """
    tlen = phy.total_branch_length()
    if target >= tlen:
        sys.stderr.write(f"prune_length: requirement already met "
                         f"({tlen}<={target}")
    tips = phy.get_terminals()

    count = 0  # Alternate by removing shorter and longer branches
    while tlen > target+range:
        # print(f"\n>>>> len: {tlen}, number of tips: {len(tips)}")
        # we have to re-sort every time because removing a branch
        # will lengthen another branch
        if count % 2 == 0 :
            tips = sorted(tips, key=lambda x: x.branch_length, reverse=True)
        else:
            tips = sorted(tips, key=lambda x: x.branch_length)
        _ = phy.prune(tips[0])
        tlen -= tips[0].branch_length
        tips = tips[1:]  # update list
        count += 1
    return phy

def run_fasttree(aln):
    # Create a tree from cleaned a nucleotide alignment
    print(f"Creating tree from: {aln}")
    p = subprocess.Popen(['fasttree', '-quiet','-nt', aln], 
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

def align_and_build_tree(clust_seqs, aln_name):
    """
    Perform a codon-aware alignment of nucleotide sequences with mafft
    Create a phylogenic tree from the alignment with fasttree
    :param clust_seqs: dict, keyed by header, 
                             values are sequences of CDSs encoding the same protein
    :aln_name: str, name of the file to store the alignment so fastree can read it
    """
    # Codon-aware sequence alignment
    before_prune_aln = align_codon_aware(clust_seqs)

    # Store alignment file so fasttree can process it
    with open(aln_name, 'w') as file:
            file.write(before_prune_aln.getvalue())
    file.close
    
    # Build tree
    tree = run_fasttree(aln_name)
    handle = StringIO(tree.decode())
    before_prune_phy = Phylo.read(handle, "newick")

    return(before_prune_aln, before_prune_phy)

def filter_seqs(grouped_seqs):
    """
    Filter sequences?
    """
    filtered = {}
    # Sequences from long branches in hepatoviruses
    temp = ["KX420952.1", "MG18943.1", "KT229612.1", "OM451167.1", 
            "KT229611.1", "OR45340.1", "OR452344.1", "OR452343.1", 
            "OQ559662.1", "KT819575.1", "EU140838.1", "D00924.1"]

    # For each protein
    for protein, seqs in grouped_seqs.items():
        print(protein)
        seq_headers = list(seqs)
        res = []
        # Loop trough all the headers looking for bad accns
        for header in seq_headers:
            match = [header for ele in temp if(ele in header)]
            if match:
                res.extend(match)
        final_headers = [header for header in seq_headers if header not in res]
        print(res)
        print(len(seq_headers))
        print(len(final_headers))
        filtered[protein] = {header: grouped_seqs[protein][header] for header in final_headers}
    
    return filtered
    
def prune_align_build_tree(clust_seqs, before_prune_phy, target, range=0.05):
    """
    Prune a phylogenetic tree to a target length, build a new codon-aware alignment and tree
    :param clust_seqs: dict, keyed by header, 
                             values are sequences of CDSs encoding the same protein
    :param before_prune_phy: Phylo object, tree built from before_pruned_aln
    :param target: str, limit for the length of your tree
    """

    target = float(target)
    # Get tree info
    # tip_names = set([tip.name for tip in phy.get_terminals()])
    tlen = before_prune_phy.total_branch_length()

    # Prune tree 
    if target < tlen:
        
        # Prune tree
        try:
            pruned_tree = prune_length(before_prune_phy, target=target, range=range)

        # Error in pruned tree
        except Exception as e:
            # It is usually not root found
            raise Exception(f"{e}")
        
        # Pruned tree is now too short
        if pruned_tree.total_branch_length() < (target - range):
            new_l = pruned_tree.total_branch_length()
            n_seqs = len(pruned_tree.get_terminals())
            raise Exception(f"Pruned tree is too short: {new_l}long, {n_seqs} seqs")

        # Get sequences after prunning
        pruned_seqs = {}
        for tip in pruned_tree.get_terminals():
            pruned_seqs[tip.name] = clust_seqs[tip.name]
        
        # Re-align sequences on tips
        pruned_aln = align_codon_aware(pruned_seqs)

    else: 
        raise Exception(f"""\n>> Target length ({target}) 
                is shorter than tree length 
                ({round(tlen, 3)})\n""")

    return (pruned_aln, pruned_tree)

def separate_clustered_seqs(sequences, clusters_info, n_prots):
    """
    Separate sequences based on their cluster asignment
    :param sequences: fasta file with CDSs
    :param clus_info: csv file with sequence header and cluster classification
    :param n_prots: number of proteins to analyse
    :return all_clus_seqs: dict, keyed by cluster, 
                                 values are dicts of sequences in cluster
                                 (header:sequence)
    """
    seqs = {rec.name: rec.seq for rec in SeqIO.parse(sequences, "fasta")}
    n_prots = args.n_prots
            
    # Read clustering results
    reader = csv.DictReader(clusters_info)
    
    # Separate sequences from clustering results
    all_clus_seqs = {} 
    for row in reader:
        cluster = row['clusters']
        label = row['name']
        clus_name = row['clus.name']
        if n_prots:
            if int(cluster) > n_prots:
                print(f" Skipping {label} from cluster: {cluster}")
                continue

        if clus_name not in all_clus_seqs:
            all_clus_seqs[clus_name] = {}
        all_clus_seqs[clus_name][label] = seqs[label]

    return all_clus_seqs

def create_report():
    report = {}

    return report

if __name__=="__main__":
    args = parse_args()

    # Sequences were clustered with kmer dist and hclust
    if args.clusters_info:
        cds_file = args.cds_file[0]
        grouped_seqs = separate_clustered_seqs(cds_file,
                                               args.clusters_info, 
                                               args.n_prots)

    # Sequences have been previously grouped
    else:
        grouped_seqs = {}
        for file in args.cds_file:
            file_name = os.path.basename(file).split('.')[0]
            grouped_seqs[file_name] = {rec.name: rec.seq for rec in SeqIO.parse(file, "fasta")}
    

    bad_trees = []  # Clusters with tree errors
    finished_analysis = []
    # Align, make tree, prune tree, measure selection
    result = {}  # Store prunning results
    for cluster in grouped_seqs:
        
        print(f"\n>>>> Processing cluster: {cluster} <<<<\n")
        clust_seqs = grouped_seqs[cluster]
        
        # Define a name for alignment file
        cluster_label = f"{args.label}_{cluster}"
        aln_name = f"{cluster_label}.codon_aln.mafft.fasta"
        
        # Get codon aware alignment with mafft and tree with fasttree
        before_prune_aln, before_prune_phy = align_and_build_tree(clust_seqs, aln_name)       
        print(f"\nStarting tree length: {before_prune_phy.total_branch_length()}\n")

        # Store results for cluster analysis
        result[cluster] = {'cluster': cluster, 
                           'initial_seqs': len(clust_seqs),
                           'pruned': False, 
                           'inital_tree_lenght': round(before_prune_phy.total_branch_length(), 3),
                           'final_seqs': 'NA',
                           'final_tree_length':'NA'}
        ks = result[cluster].keys()
        
        # Save phylogeny and codon aware alignment before prunning for debugging
        if args.save_before:
            with open(f"{cluster_label}_before_prun.fasta", 'w') as file:
                file.write(before_prune_aln.getvalue())
            file.close()
            Phylo.write(before_prune_phy, f"{cluster_label}.before_prun.tree", 'newick')

        # Prune tree
        if args.prune:
            try:
                pruned_aln, pruned_tree = prune_align_build_tree(clust_seqs,
                                                before_prune_phy, 
                                                args.prune)
                
                # Save phylogeny
                Phylo.write(pruned_tree, f"{cluster_label}.tree", 'newick')
                # Overwrite alignment with pruned_tree sequences
                with open(aln_name, 'w') as file:
                    file.write(pruned_aln.getvalue())
                
                # Update results when prune successfully
                result[cluster]['pruned'] = True
                result[cluster]['final_seqs'] = len(pruned_tree.get_terminals())
                result[cluster]['final_tree_length'] = round(pruned_tree.total_branch_length(), 3)

                finished_analysis.append(cluster)
                
            # Error in pruned tree
            except Exception as e:
                print(f"\n-----------------------------------------------")
                print(f"\tError while pruning tree of cluster {cluster}:")
                print(f"\t'{e}'")
                print(f"\tConsider removing sequences from longest branches")
                print(f"-------------------------------------------------\n") 
                bad_trees.append(cluster)
                continue

        # Measure selection with FUBAR
        if args.run_sel:
            run_selection_pipeline(aln_name,
                                   cluster_label)
    
    # Save prunning results in report file 
    if args.save_prune_res:
        header = list(result[list(result.keys())[0]].keys())
        with open(f'{args.label}.report.csv', 'w') as csvfile:
            fieldnames = header
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            
            writer.writeheader()
            for clust in result.values():
                writer.writerow(clust)

    print(f"\n>>> Unsuccessful analysis for clusters: {bad_trees}")
    print(f">>> Successful analysis for clusters: {finished_analysis}\n")