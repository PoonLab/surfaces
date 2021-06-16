#find_rev.py reverse complement of nucleotide anti-sense seqeucnes 
#if (-1) in the header then reverse complement of that fasta entry is written out to a fasta file in the rev_directory
#input directory to fasta file

import os
import argparse
from glob import glob
import subprocess

directory_arg = '/home/sareh/data/ref_cds'
rev_directory_arg = '/home/sareh/data/rev_nuc_refseq_CDS'

complement_dict = {'A':'T', 'C':'G', 'G':'C', 'T':'A',
                    'W':'S', 'R':'Y', 'K':'M', 'Y':'R', 'S':'W', 'M':'K',
                    'B':'V', 'D':'H', 'H':'D', 'V':'B',
                    '*':'*', 'N':'N', '-':'-'}

def reverse_and_complement(seq):
    rseq = seq[::-1]
    rcseq = ''
    for i in rseq:  # reverse order
        rcseq += complement_dict[i]
    return rcseq


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

def main():
    count = 0
    os.mkdir(rev_directory_arg)
    for dir in os.listdir(directory_arg):
        new_dir_path=os.path.join(rev_directory_arg, dir)
        os.mkdir(new_dir_path)
        #check empty directories and erase them
        dir_path=os.path.join(directory_arg,dir)
        for file_name in (os.listdir(dir_path)):
            file_path =os.path.join(dir_path,file_name)
            with open(file_path, 'r') as f:
                for h, sequence in iter_fasta(f):
                    term =(h.split(","))[3]
                    if term == "-1":
                        with open(os.path.join(new_dir_path, file_name),'w') as outfile:
                            revseq=reverse_and_complement(sequence)
                            outfile.write(">{}\n{}\n".format(h,revseq))

if __name__ == '__main__':
    main()
