"""
Combine testminimap and cutter_mpi.py
mkdir cut_cds | mkdir cutter_scores | mkdir gapped
"""

from glob import glob
import tempfile
import subprocess
import re
import os
import sys
from datetime import datetime

# mpi
from mpi4py import MPI
the_world = MPI.COMM_WORLD
my_number = the_world.Get_rank()
total_number = the_world.Get_size()

def get_accn(file):
    """
    takes in file path
    Stream through txt file retrieve accns into a list
    """
    lst_accns = []
    with open(file, 'r') as f:
        for line in f:
            line = line.strip('\n')
            lst_accns.append(line)
    a = tuple(lst_accns)
    return(a)

def convert_fasta (handle):
    """
    takes in handel to an open fasta file 
    """
    result = []
    sequence = ''
    # handle = open(file,'r') 
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

def pdist(s1, s2):
    """
    Calculate the p-distance (proportion of nucleotide differences between
    two sequences)
    """
    ndiff = 0
    for i, nt1 in enumerate(s1):
        if nt1 in 'N-':
            continue
        nt2 = s2[i]
        if nt1 != nt2 and nt2 not in 'N-':
            ndiff += 1
    return ndiff/len(s1)

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
    

def minimap2(query, refseq, nthread=4, mm2bin='minimap2'):
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

def mafft(query, ref, trim=True):
    handle = tempfile.NamedTemporaryFile(delete=False)
    s = '>ref\n{}\n>query\n{}\n'.format(ref, query)
    handle.write(s.encode('utf-8'))
    handle.close()
    
    # call MAFFT on temporary file
    stdout = subprocess.check_output(['mafft', '--quiet', handle.name]) #Failed to align Error
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


def cutter_minimap(ref, fasta, outfile, csvfile, gapped):
    """
    :param ref:  open stream in read mode to FASTA file containing reference
                 gene sequence
    :param fasta:  open stream in read mode to FASTA file containing query
                   genome sequences
    :param outfile:  open stream in write mode to output trimmed sequences
                     in FASTA format
    :param csvfile:  open stream in write mode to output alignment scores
    
    :param gapped:  open stream in write mode to output ref seq that mafft inserted gaps into
    """
  
    csvfile.write('header,align.score\n')
    gapped.write("header,rgene_mafft,qgene_mafft,qgene_nogap,gap_qgene_mafft,gap_rgene_mafft,original_rgene\n")    
    
    # parse query genomes from input FASTA
    queries = convert_fasta(fasta)
    
    # the reference cds seqeunce
    original_rgene = convert_fasta(ref)[0][1]
    
    for qh, qs in queries:
        # minimap2 alignment: use minimap2 to extract homologous gene in query genome
        qgene_minimap = minimap2(query=qs, refseq=original_rgene)

        if qgene_minimap is None:
            # failed to align reference gene to query genome - no homology?
            sys.stdout.write("{}\n".format(qh))
            # outfile.write('>{}\n{}\n'.format(qh, qgene_minimap))
            continue

        # pairwise alignment
        qgene_mafft, rgene_mafft = mafft(query=qgene_minimap, ref=original_rgene, trim=False)
        
        # check for gaps in reference gene (insertion in query)       
               # qgene_nogap = re.sub('-','',qgene_mafft)

        qgene_nogap = ""
        for i, rn in enumerate(rgene_mafft):
            if rn == '-':
                # skip insertion relative to reference
                continue
            qgene_nogap += qgene_mafft[i] 

        outfile.write(">{}\n{}\n".format(qh, qgene_nogap))

        # Print length of each gene in the gapped outfile 
        gap_qgene_mafft = qgene_mafft.count("-")
        gap_rgene_mafft = rgene_mafft.count("-")
        gapped.write("{},{},{},{},{},{},{}\n".format(qh,str(len(rgene_mafft)),str(len(qgene_mafft)),str(len(qgene_nogap)),str(gap_qgene_mafft),str(gap_rgene_mafft),len(original_rgene)))
#        print("rgene_mafft" + str(len(rgene_mafft)))
#        print("qgene_mafft" + str(len(qgene_mafft)))
#        print("qgene_nogap"  + str(len(qgene_nogap))) 
#        print("gap_qgene_mafft" + str(gap_qgene_mafft))
#        print("gap_rgene_mafft" + str(gap_rgene_mafft))str(gap_qgene_mafft)

        # calculate p-distance as a measure of similarity/alignment score
        p = pdist(qgene_mafft, rgene_mafft)
        csvfile.write('{},{:1.3f}\n'.format(qh, 100*p))

def test():
    count = 0

    ref_gene_directory = "/home/sareh/hiv/NC_001802"
    outdir = "/home/sareh/hiv/cds/100"

    # Path to query fasta file with all the full genomes 
    fasta = "/home/sareh/hiv/prune/100_pruned/pruned_NC_001802.fa"

    for file in os.listdir(ref_gene_directory):
            count += 1
            if count % total_number != my_number:
                continue

            # Path to FASTA with gene sequence from Reference genome
            ref = os.path.join(ref_gene_directory, file)

            # Path to write FASTA output (cut gene sequences)
            outfile = ("{}/cut_cds/cutter_cds_{}".format(outdir, file))

            # Path to write CSV output (alignment scores)
            csvfile = ("{}/cutter_scores/cutter_scores_{}".format(outdir, file))

            #path to write out rgenes where gap was inserted
            gapped = ("{}/cds_gapped/{}".format(outdir,file))

            # run analysis
            outfile_handle = open(outfile, 'w+')  # append
            csvfile_handle = open(csvfile, 'w+')
            ref_handle = open(ref)
            fasta_handle = open(fasta)
            gapped_handle = open(gapped, 'w+')

            cutter_minimap(ref_handle, fasta_handle, outfile_handle, csvfile_handle,gapped_handle)

            ref_handle.close()
            fasta_handle.close()
            outfile_handle.close()
            csvfile_handle.close()
            gapped_handle.close()

def main():
    count = 0
    
    # Directories 
    outdir = '/home/sareh/data/sequences/cut_cds' # where all the files would be written into
    ref_home_directory = '/home/sareh/surfaces/find_cds/corrected_ref_cds'
    query_directory = '/home/sareh/data/sequences/pruned_genome'
    accession_file = '/home/sareh/surfaces/find_cds/virus_pruned_genome.txt'

    accns = get_accn(accession_file) # all of the files to examine

    for accn in accns:
        ref_gene_directory = os.path.join(ref_home_directory, accn)
        
        # Path to FASTA with genome sequences to process
        #fasta = ('{}/Pruned_nuc_{}'.format(query_directory, accn))
        fasta = ('{}/Pruned_nuc_{}'.format(query_directory, accn))
        
        for file in os.listdir(ref_gene_directory):
            count += 1
            if count % total_number != my_number:
                continue

            # Path to FASTA with gene sequence from Reference genome
            ref = os.path.join(ref_gene_directory, file)

            # Path to write FASTA output (cut gene sequences)
            outfile = ("{}/cut_cds/cutter_cds_{}".format(outdir, file))

            # Path to write CSV output (alignment scores)
            csvfile = ("{}/cutter_scores/cutter_scores_{}".format(outdir, file))

            #path to write out rgenes where gap was inserted
            gapped = ("{}/gapped/{}".format(outdir,file))

            #os.path.isfile(path) | os.path.exists(path) | path.exists
            if os.path.isfile(outfile):
                # output file exists, skip to next job
                continue

            # progress monitoring
            sys.stdout.write("[{} {}/{}] starting job {} {}\n".format(
                datetime.now().isoformat(), my_number, total_number, count, file))
            sys.stdout.flush()

            # run analysis
            outfile_handle = open(outfile, 'w+')  # append
            csvfile_handle = open(csvfile, 'w+')
            ref_handle = open(ref)
            fasta_handle = open(fasta)
            gapped_handle = open(gapped, 'w+')

            cutter_minimap(ref_handle, fasta_handle, outfile_handle, csvfile_handle,gapped_handle)

            #if os.path.isfile(outfile) and count % total_number == my_number:
            ref_handle.close()
            fasta_handle.close()
            outfile_handle.close()
            csvfile_handle.close()
            gapped_handle.close()

if __name__ == "__main__":
    test()
