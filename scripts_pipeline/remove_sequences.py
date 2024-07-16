description = """
Remove any sequences from target FASTA that are not in reference FASTA.
This is used after manually removing problematic sequences from a 
FASTA file, and when there is a second file with the same headers.
"""

from Bio import SeqIO
import sys
import argparse


parser = argparse.ArgumentParser(description=description)
parser.add_argument("reffile", type=argparse.FileType('r'), 
                    help="Reference FASTA file with headers manually removed.")
parser.add_argument("target", type=argparse.FileType('r'),
                    help="Target FASTA file for modification.")
parser.add_argument("-o", "--outfile", type=argparse.FileType('w'), default=sys.stdout,
                    help="Output file, defaults to stdout.")
args = parser.parse_args()

reffile = SeqIO.parse(args.reffile, "fasta")
headers = set([record.description for record in reffile])
missing = []
for record in SeqIO.parse(args.target, "fasta"):
    if record.description in headers:
        args.outfile.write(record.format("fasta"))
        continue
    missing.append(record.description)

if missing:
    sys.stderr.write(f"Skipped {len(missing)} records from target file.\n")

