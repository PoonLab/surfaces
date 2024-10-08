from Bio import SeqIO
import sys
import argparse
import re

description = """
Split a FASTA file containing two sets of non-homologous gene sequences.
For instance, one of the segments of the influenza A virus genome encodes
M1 and M2 genes, where M2 is produced by splicing a 5' portion of M1 with 
downstream sequence.  It can be differentiated from M1 by the location 
field.
"""
pat = re.compile("[<>]?([0-9]+)_[<>]?([0-9]+)")

parser = argparse.ArgumentParser(description = "Split FASTA containing M1 and M2 genes "
                                 "into separate FASTA files")
parser.add_argument("input", type=argparse.FileType('r'),
                    help="input FASTA file")
parser.add_argument("out1", type=argparse.FileType('w'),
                    help="path to write FASTA of first sequences ('.' in loc / >maxlen)")
parser.add_argument("out2", type=argparse.FileType('w'),
                    help="path to write FASTA of second sequences (<=maxlen)")
parser.add_argument("--maxlen", type=int, default=None,
                    help="optional, partition sequences by length")
args = parser.parse_args()

for record in SeqIO.parse(args.input, 'fasta'):
    loc = record.description.split('-')[-1]
    check = '.' in loc
    if args.maxlen:
        check = len(record.seq) > args.maxlen
    
    SeqIO.write(record, args.out1 if check else args.out2, 'fasta')

