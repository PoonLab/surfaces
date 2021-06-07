#The reverse complement of nucleotide sequences in a fasta file

import os
import argparse
from Bio.Seq import Seq

def iter_fasta(handle):
    """
    Parse open file as FASTA.  Returns a generator
    of handle, sequence tuples.
    :param handle:  open stream to FASTA file in read mode
    :yield tuples, (header, sequence)
    """
    h, sequence = None, ''
    for line in handle:
        if line.startswith('>'):
            if len(sequence) > 0:
                yield h, sequence
                sequence = ''
            h = line.lstrip('>').rstrip()
        else:
            sequence += line.strip().upper()
    yield h, sequence

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--file', type=str,
                       help='path to fasta file with nucleotide sequences')
    parser.add_argument('--outdir', default='/home/sareh/data', type=str,
                        help='path to directory to write outputs.')
    return parser.parse_args()

def main():
    args = parse_args()
    with open (args.file, 'r') as file:
        out_path = os.path.join(args.outdir, os.path.basename(args.file))
        with open(out_path, 'w') as outfile:
            for h,sequence in iter_fasta(file):
                seq=Seq(sequence)
                new_seq = seq.reverse_complement()
                outfile.write(">{}\n{}\n".format(h,new_seq))

if __name__ == '__main__':
    main()
