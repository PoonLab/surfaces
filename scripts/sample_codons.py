description = """
Sample codon sites from an alignment at random with or without replacement.
Write the result to a new file.
"""

from Bio import AlignIO
import argparse
import random


parser = argparse.ArgumentParser(description=description)
parser.add_argument("infile", type=argparse.FileType('r'), help="input, path to FASTA file")
parser.add_argument("--replace", '-r', action="store_true", 
                    help="option, set to sample with replacement")
parser.add_argument("-f", "--format", type=str, default='fasta',
                    help="option, sequence file format (default: fasta)")
parser.add_argument("-n", type=int, default=1, 
                    help="option, number of replicates (default: 1)")
parser.add_argument("-p", "--prefix", type=str, 
                    help="prefix (no '.fasta' extension) for path to write output."
                    " If -n is greater than 1, '_n' will be appended to prefix.")
args = parser.parse_args()

assert args.n > 0, "Error, -n must be greater than zero!"

records = AlignIO.read(args.infile, args.format)
