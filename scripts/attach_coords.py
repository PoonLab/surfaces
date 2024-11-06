from Bio import SeqIO
import csv
import os
from glob import glob
from consensus import get_columns, conseq
import tempfile
import sys


# parse annotation CSV
refseq = {}
handle = open("data/annotation.csv")
reader = csv.DictReader(handle)
for row in reader:
    virus = row['virus']
    if virus not in refseq:
        refseq.update({virus: {}})
    refseq[virus].update({row['protein']: row['refseq']})


# generate step3 consensus sequences
files = glob("data/3_codaln/*_step3.fasta")
for fn in files:
    tokens = os.path.basename(fn).split('_')
    virus = tokens[0]
    if virus not in refseq:
        sys.stderr.write(f"ERROR: Could not find virus {virus} in annotations\n")
        sys.exit()
    
    protein = ' '.join(tokens[1:-1])
    if protein not in refseq[virus]:
        sys.stderr.write(f"ERROR: Could not find {virus} protein '{protein}' in annotations\n")
        sys.exit()
    
    aln = SeqIO.parse(fn, 'fasta')
    columns = get_columns(aln, is_nuc=True)
    cseq = conseq(columns, thresh=0.5, is_nuc=True)
    break