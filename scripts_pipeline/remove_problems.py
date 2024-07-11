"""
Remove problematic sequences from a FASTA file that are listed in 
another plain text file
"""

from Bio import SeqIO
import sys


# load sequences from FASTA
handle = open("seg5_CDSs.fasta") 
records = SeqIO.parse(handle, 'fasta')

seqdict = dict([(record.description, record) for record in records])

# load list of problematic sequences
handle = open("seg5_problem.txt")
for line in handle:
    label = line.strip()
    if len(label) == 0:
        continue

    found = False
    for header, record in seqdict.items():
        if header.startswith(label):
            print(f"Found match for {label}")
            found = True
            break

    if found:
        seqdict.pop(header)
    else:
        print(f"label {label} not found")
        sys.exit()


# export remaining sequences to new file
outfile = open("no_problem.fasta", 'w')
for header, record in seqdict.items():
    outfile.write(f">{header}\n{record.seq}\n")

outfile.close()

