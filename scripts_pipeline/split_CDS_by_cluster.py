import csv
import argparse
import os

from Bio import SeqIO

description = """\
Create separate fasta file based on cluster classification.
"""

def parse_args():
    parser = argparse.ArgumentParser(
        description="Create multiple fasta, one for each cluster of proteins"
    )
    parser.add_argument(
        'cds_file', type=str,
        help = 'Path to fasta file with CDSs'
    )
    parser.add_argument(
        '--clusters_info', '-ci', type=argparse.FileType('r'), default=False,
        help='Path to file containing cluster output in csv format'
    )
    parser.add_argument(
        '--location', '-l' , type=str, default='cluster',
        help = 'location for the CDS files'
    )
    parser.add_argument(
        '-n', '--n_prots', type=int, default=False, 
        help = 'Maximum number of clusters allowed'
    )

    return parser.parse_args()

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


if __name__=="__main__":
    args = parse_args()

    # Sequences were clustered with kmer dist and hclust
    if args.clusters_info:
        cds_file = args.cds_file
        dir = os.path.dirname(args.cds_file)
        grouped_seqs = separate_clustered_seqs(cds_file,
                                               args.clusters_info, 
                                               args.n_prots)

        for clus_name, seqs in grouped_seqs.items():
            with open(f"{args.location}/{clus_name}.fasta", 'w') as file:
                for desc, seq in seqs.items():
                    file.write(f">{str(desc)}\n{str(seq)}\n")
