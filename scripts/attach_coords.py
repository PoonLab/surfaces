from Bio import SeqIO
import csv
import os
from glob import glob
from Bioplus.consensus import conseq
from Bioplus.pair_align import mafft
import tempfile
import sys


# parse annotation CSV
accns = {}
handle = open("data/annotation.csv")
reader = csv.DictReader(handle)
for row in reader:
    virus = row['virus']
    if virus not in accns:
        accns.update({virus: {}})
    accn = row['refseq'].split('.')[0]
    accns[virus].update({row['protein']: accn})


# retrieve sequences from Genbank files
refseqs = {}
for fn in glob("data/refseq/*.gb"):
    record = SeqIO.read(fn, 'genbank')
    accn = os.path.basename(fn).split('.')[0]
    refseqs.update({accn: record.seq})


# process step 3 files
files = glob("data/3_codaln/*_step3.fasta")
for fn in files:
    # map to accession
    tokens = os.path.basename(fn).split('_')
    virus = tokens[0]
    if virus not in accns:
        sys.stderr.write(f"ERROR: Could not find virus {virus} in annotations\n")
        sys.exit()
    
    protein = ' '.join(tokens[1:-1])
    if protein not in accns[virus]:
        sys.stderr.write(
            f"ERROR: Could not find {virus} protein '{protein}' in annotations\n")
        sys.exit()

    accn = accns[virus][protein]

    # retrieve reference
    refseq = refseqs.get(accn, None)
    if refseq is None:
        sys.stderr.write(f"WARNING: failed to retrieve reference for {accn}\n")
        continue

    # generate consensus sequence
    aln = SeqIO.parse(fn, 'fasta')
    cseq = conseq(aln, thresh=0.5, is_nuc=True)
      
    # align to corresponding reference
    aquery, aref = mafft(cseq, refseq)
    
