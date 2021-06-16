import subprocess
import tempfile
import argparse
#from seqUtils import get_boundaries, convert_fasta
import sys
import re
import os
from gotoh2 import Aligner

aligner = Aligner()  # default settings


def convert_fasta (handle):
    result = []
    sequence = ''
    for line in handle:
        if line.startswith('$'): # skip header line
            continue
        elif line.startswith('>') or line.startswith('#'):
            if len(sequence) > 0:
                result.append([h,sequence])
                sequence = ''   # reset
            h = line.strip('>#\n')
        else:
            sequence += line.strip('\n').upper()
            
    result.append([h,sequence]) # handle last entry
    return result


def get_boundaries (str):
    gap_prefix = re.compile('^[-]+')
    gap_suffix = re.compile('[-]+$')
    # return a tuple giving indices of subsequence without gap prefix and suffix
    res = [0,len(str)]
    left = gap_prefix.findall(str)
    right = gap_suffix.findall(str)
    if left:
        res[0] = len(left[0])

    if right:
        res[1] = len(str) - len(right[0])

    return res

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('ref', type=argparse.FileType('r'), 
                        help="FASTA file with reference sequence to cut out of queries")
    parser.add_argument('fasta', type=argparse.FileType('r'), 
                        help="FASTA file containing query sequences")
    parser.add_argument('outfile', type=argparse.FileType('w'), default=None,
                        help="path to write output")
    parser.add_argument('--echo', action='store_true', help="report alignment scores")
    return parser.parse_args()


def mafft(query, ref):
    handle = tempfile.NamedTemporaryFile(delete=False)
    s = '>ref\n{}\n>query\n{}\n'.format(ref, query)
    handle.write(s.encode('utf-8'))
    handle.close()
    
    # call MAFFT on temporary file
    stdout = subprocess.check_output(['mafft', '--quiet', handle.name])
    stdout = stdout.decode('utf-8')
    result = convert_fasta(stdout.split('\n'))
    aligned_ref = result[0][1]
    aligned_query = result[1][1]
    
    # trim aligned query sequence to extent of reference
    left, right = get_boundaries(aligned_ref)
    trimmed_query = aligned_query[left:right]
    os.remove(handle.name)  # clean up
    return(trimmed_query)


def gotoh(query, ref, echo=False):
    aref, aquery, ascore = aligner.align(ref, query)
    left, right = get_boundaries(aref)
    if echo:
        print(ascore / (right-left))
        trimmed_query = aquery[left:right]
        return trimmed_query.replace('-', '')  # remove all gaps

def cutter():
    args = parse_args()
    
    # read reference sequence from file
    refseq = convert_fasta(args.ref)[0][1]

    fasta = convert_fasta(args.fasta)
    for h, s in fasta:
        query = s.replace('-', '')
        #trimmed = mafft(query, refseq)
        genome_ID = h.split(",")[1] 
        print("{}".format(genome_ID)) 
        #print header with virusID and genomeID
        #print("{}".format(h)) 
        trimmed = gotoh(query, refseq, echo=args.echo)
        args.outfile.write('>{}\n{}\n'.format(h, trimmed))

if __name__ == '__main__':
    cutter()
