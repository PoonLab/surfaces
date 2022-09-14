import os
import sys
import argparse

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--indir', type=str,
                       help='directory of tmhmm 3 line output')
    parser.add_argument('--outfile', type=str,
                        help='tmhmm csv output')
    return parser.parse_args()


def parse_tmhmm_output(handle):
    """
    Parse open file as FASTA.  Returns a generator
    of handle, sequence tuples.
    :param handle:  open stream to FASTA file in read mode
    :yield tuples, (header, sequence)
    """
    h, sequence = None, ''
    for num, line in enumerate(handle):
        # first line is the header 
        if num == 1:
            h = line.strip(">")
        # secound line is the aa seq 
        elif num == 2:
            seq = line.strip("\n")
        # third line is the predictions 
        elif num == 3" 
            prdct = line.strip("\n")
            
    return(h, seq, prdct)


def main():
    
    args = parse_args()
    for accn in os.listdir(args.indir):
        infile = os.path.join(args.indir, accn)
        in_handle = open(infile)
        h, seq, prdct = parse_tmhmm_output(handle)
        print(h)
        print(seq)
        print(prdct)
        
if __name__ == '__main__':
    main()
