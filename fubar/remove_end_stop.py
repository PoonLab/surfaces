#!/usr/bin/python

import re
import os
import sys

file_path = sys.argv[1] #fasta file
out_path = sys.argv[2] #fasta file

start = ["ATG"] 
stop = ["TAA","TAG","TGA"]

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


def remove_stop(seq):
    """
    takes in a string of sequences and removes stop codon only from the end
    """
    last_codon = seq[-3:]
    if last_codon in stop:
        new_seq = seq[:-3]
    else: 
        new_seq = seq
    return(new_seq)

def main():
    file_handle = open(file_path,'r')
    out_handle = open(out_path,'w')

    for h,s in iter_fasta(file_handle):
        new_seq = remove_stop(s)
        out_handle.write(">{}\n{}\n".format(h,new_seq))

    file_handle.close()
    out_handle.close()

if __name__ == '__main__':
    main()
