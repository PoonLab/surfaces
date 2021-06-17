#cutter mpi
#python3 cutter.py /home/sareh/surfaces/retrieve_sequences/cutter_test/cutter_ref_example.txt
# /home/sareh/data/Pruned_sequence/Pruned_nuc_NC_001547
# /home/sareh/surfaces/retrieve_sequences/cutter_test/output_NC_00154 --echo
# > /home/sareh/surfaces/retrieve_sequences/cutter_test/NC_001547_score.txt

from mpi4py import MPI
from gotoh2 import Aligner
import argparse
from glob import glob
import subprocess
import tempfile
import os
import re
import sys

the_world = MPI.COMM_WORLD
my_number = the_world.Get_rank()
total_number = the_world.Get_size()

directory = '/home/sareh/data/ref_cds'

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


def gotoh(query, ref):
    aref, aquery, ascore = aligner.align(ref, query)
    left, right = get_boundaries(aref)
    ascore = ascore / (right-left)

    trimmed_query = aquery[left:right].replace('-', '')  # remove all gaps
    return trimmed_query, ascore


def cutter(ref, fasta, outfile, csvfile):
    """
    :param ref:  open stream in read mode to FASTA file containing reference
                 gene sequence
    :param fasta:  open stream in read mode to FASTA file containing query
                   genome sequences
    :param outfile:  open stream in write mode to output trimmed sequences
                     in FASTA format
    :param csvfile:  open stream in write mode to output alignment scores
    """
    # read reference sequence from file
    refseq = convert_fasta(ref)[0][1]

    csvfile.write('header,align.score\n')

    fasta = convert_fasta(fasta)
    for h, s in fasta:
        # remove gaps in query before attempting alignment
        query = s.replace('-', '')

        #trimmed = mafft(query, refseq)
        genome_ID = h.split(",")[1] 
        print("{}".format(genome_ID)) 

        trimmed, ascore = gotoh(query, refseq)

        outfile.write('>{}\n{}\n'.format(h, trimmed))
        csvfile.write('"{}",{}\n'.format(h, ascore))


if __name__ == '__main__':
    count = 0

    for dir in os.listdir(directory):
        if 'NC_' in dir:
            #path = ("{}/{}".format(directory,dir)) #print(dir) #NC_044047
            path = os.path.join(directory, dir)
            for file in os.listdir(path):
                count += 1 
                if count % total_number == my_number:
                    #print(file) #NC_044047_YP_009679042.1
                    ref = ("{}/{}/{}".format(directory, dir, file))
                    fasta = ("/home/sareh/data/Pruned_sequence/Pruned_nuc_{}".format(dir))
                    outfile = ("/home/sareh/cutter/rev_cutter/cutter_cds/rev_cutter_{}".format(file))
                    score = ("/home/sareh/cutter/rev_cutter/rev_cutter_scores/{}".format(file))

                    cutter(ref, fasta, outfile, score) 
                    break

