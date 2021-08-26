#! /usr/bin/python3

"""
retrieve sequences labelled as "complete genome"
"""
import argparse

#fasta file from NCBI
input_path = args[1]
#complete_fasta
output_path = args[2]

def iter_fasta (handle):
    """
    Parse open file as FASTA.  Returns a generator
    of handle, sequence tuples.
    """
    sequence = ''
    for i in handle:
        if i[0] == '$': # skip h info
            continue
        elif i[0] == '>' or i[0] == '#':
            if len(sequence) > 0:
                yield h, sequence
                sequence = ''   # reset containers
            h = i.strip('\n')[1:]
        else:
            sequence += i.strip('\n').upper()
    yield h, sequence

def main():
    input_handle = open(input_path,'r')
    outfile = open(output_path,"w")
    count = 0
    for h,s in iter_fasta(input_handle):
        if "complete genome" in h and "nearly" not in h:
            outfile.write(">{}\n{}\n".format(h,s))
            count += 1 

    print(count)
    input_handle.close()
    outfile.close()

if __name__ == '__main__':
    main()
