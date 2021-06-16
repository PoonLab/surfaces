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
    not_zero = 0 
    not_triple = 0 #number of sequence entries in fasta file not divisible by 3 
    not_lst = [] #list of headers of fasta entries in each file not divisible by 3 
    with open(args.file, 'r') as file:
        file_name = args.file[35:]
        for h,sequence in iter_fasta(file):
            divisible = len(sequence) % 3
            count += 1 
            if divisible != 0:
                not_lst.append(h)
                not_triple += 1
        if not_triple != 0:
            not_zero+=1
        #HEARED:not_tripple,total,file_name,length_not_tripple,remainder
        print("{},{},{},{}".format(not_triple,count,file_name,len(not_lst),divisible))
    print(not_zero)
if __name__ == '__main__':
    main()
