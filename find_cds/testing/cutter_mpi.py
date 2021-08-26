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
from datetime import datetime
import subprocess


the_world = MPI.COMM_WORLD
my_number = the_world.Get_rank()
total_number = the_world.Get_size()

aligner = Aligner()  # default settings


def convert_fasta (handle):
    #takes in opened file handle
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
            
    result.append([h, sequence])  # handle last entry
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


def mafft(query, ref, trim=True):
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
    if trim:
        left, right = get_boundaries(aligned_ref)
        aligned_query = aligned_query[left:right]
    os.remove(handle.name)  # clean up
    return(aligned_query, aligned_ref)


def apply_cigar(seq, cigar):
    """
    Use CIGAR to pad sequence with gaps as required to
    align to reference.  Adapted from http://github.com/cfe-lab/MiCall
    """
    is_valid = re.match(r'^((\d+)([MIDNSHPX=]))*$', cigar)

    if not is_valid:
        raise RuntimeError('Invalid CIGAR string: {!r}.'.format(cigar))
    tokens = re.findall(r'  (\d+)([MIDNSHPX=])', cigar, re.VERBOSE)
    aligned = ''
    left = 0
    for length, operation in tokens:
        length = int(length)
        if operation in 'M=X':
            aligned += seq[left:(left+length)]
            left += length
        elif operation == 'D':
            aligned += '-'*length
        elif operation in 'SI':
            left += length  # soft clip
    return aligned


def minimap2(query, refseq, nthread=2, mm2bin='minimap2'):
    """
    :param query:  str, query genome sequence
    :param refseq:  str, reference gene sequence
    :param nthread:  int, number of threads to run minimap2
    :param mm2bin:  str, path to minimap2 binary
    """
    # create temporary FASTA file with query genome
    handle = tempfile.NamedTemporaryFile(delete=False)
    s = '>query_as_ref\n{}\n'.format(query)
    handle.write(s.encode('utf-8'))
    handle.close()

    # call minimap2
    cmd = [mm2bin, '-t', str(nthread), '-a', '--eqx', handle.name, '-']
    p = subprocess.Popen(cmd, encoding='utf8',
                         stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    output, _ = p.communicate('>test\n{}\n'.format(refseq))
    for line in output.split('\n'):
        if line == '' or line.startswith('@'):
            # split on \n leaves empty line; @ prefix header lines
            continue
        qname, flag, rname, rpos, _, cigar, _, _, _, seq = \
            line.strip('\n').split('\t')[:10]

        if rname == '*' or ((int(flag) & 0x800) != 0):
            # did not map, or supplementary alignment
            print("Failed to align")
            continue

        # validate CIGAR string
        is_valid = re.match(r'^((\d+)([MIDNSHPX=]))*$', cigar)
        if not is_valid:
            raise RuntimeError('Invalid CIGAR string: {!r}.'.format(cigar))

        rpos = int(rpos) - 1  # convert to 0-index
        aligned = apply_cigar(seq, cigar)  # check if reference gene has insertions
        return query[rpos:(rpos+len(aligned))]


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

    out_file = open(outfile, 'w+')  # append
    scorefile = open(csvfile, 'w+')

    scorefile.write('header,align.score\n')

    
    fasta = convert_fasta(fasta)
    for h, s in fasta:
        # remove gaps in query before attempting alignment
        query = s.replace('-', '')

        #trimmed = mafft(query, refseq)
        genome_ID = h.split(",")[1] 
        # print("{}".format(genome_ID)) 
        sys.stdout.write("[{} {}/{}] aligning {}\n".format(
            datetime.now().isoformat(), my_number, total_number, genome_ID))
        sys.stdout.flush()

        trimmed, ascore = gotoh(query, refseq)

        scorefile.write('{}\,{}\n'.format(h, ascore))
        out_file.write('>{}\n{}\n'.format(h, trimmed))

    out_file.close()
    scorefile.close()


def main():
    count = 0
    dir_count = 0
    file_count = 0
    outdir = '/home/sareh/surfaces/find_cds'

    # directory = '/home/sareh/data/ref_cds'
    directory = '/home/sareh/surfaces/find_cds/corrected_ref_cds'

    for dir in os.listdir(directory):
        if 'NC_' not in dir:
            continue

        dir_count += 1
        #print("dir_count is {}".format(dir_count))
        #path = ("{}/{}".format(directory,dir)) #print(dir) #NC_044047
        path = os.path.join(directory, dir)
        for file in os.listdir(path):
            count += 1
            if count % total_number != my_number:
                continue

            # path to FASTA with gene sequence from reference genome
            ref = ("{}/{}/{}".format(directory, dir, file))

            # path to FASTA with genome sequences to process
            fasta = ("/home/sareh/data/pruned_genome/Pruned_nuc_{}".format(dir))

            # path to write FASTA output (cut gene sequences)
            outfile = ("{}/cut_cds/cutter_cds_{}".format(outdir, file))

            # path to write CSV output (alignment scores)
            csvfile = ("{}/cutter_scores/cutter_scores_{}".format(outdir, file))

            #os.path.isfile(path) | os.path.exists(path) | path.exists
            if os.path.isfile(outfile):
                # output file exists, skip to next job
                continue

            # progress monitoring
            sys.stdout.write("[{} {}/{}] starting job {} {}\n".format(
                datetime.now().isoformat(), my_number, total_number, count, file))
            sys.stdout.flush()

            # run analysis
            ref_handle = open(ref)
            fasta_handle = open(fasta)
            cutter(ref_handle, fasta_handle, outfile, csvfile)
            #if os.path.isfile(outfile) and count % total_number == my_number:
            ref_handle.close()
            fasta_handle.close()
        #break


if __name__ == '__main__':
    main()

