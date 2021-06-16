#! /usr/bin/python3
#reads in fasta file and standard out headers of fasta enties where the seqeunces is not divisible by 3

import argparse

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
    return parser.parse_args()

def main ():
    args = parse_args()
    count = 0
    not_triple = 0
    with open(args.file, 'r') as file:
        print(args.file)
        for h,sequence in iter_fasta(file):
            divisible = len(sequence) % 3
            count += 1 
            if divisible != 0:
                print("{}|{}".format(divisible,h))
                not_triple += 1
        print("Total:{}".format(count))
        print("Not:{}".format(not_triple))

if __name__ == '__main__':
    main()
