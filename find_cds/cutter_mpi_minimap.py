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
from operator import itemgetter
from itertools import groupby
from pathlib import Path

# all dir
outdir = "/home/sareh/93/cut_cds_pruned_nbr_genome"
ref_home_directory ="/home/sareh/93/ref_cds/"
query_directory = "/home/sareh/93/pruned_nbr_genome/"
accession_file ="/home/sareh/93/ref_accns93.txt"

"""
# ncbi directories
outdir = "/home/sareh/2020/july/cut_cds/"
ref_home_directory ="/home/sareh/2020/july/ref_cut_cds/"
query_directory = "/home/sareh/2020/july/pruned_nbr_genome/"
accession_file ="/home/sareh/2020/july/pruned_ref_accn.txt"

# lin Directories
outdir = "/home/sareh/2020/sequences/lin/lin_cut_cds"  # where all the files would be written into
ref_home_directory = "/home/sareh/2020/sequences/lin/lin_ref_cds/" #the ref_cds
query_directory = '/home/sareh/2020/sequences/lin/lin_nuc_genome/'
#accession_file = '/home/sareh/2020/comp/test/lin_test.txt'
accession_file = '/home/sareh/2020/sequences/lin/lin_ref_accn.txt'

# TEST Directories
outdir = "/home/sareh/2020/sequences/ncbi/test/test_json/"
ref_home_directory ="/home/sareh/2020/sequences/ncbi/ncbi_ref_cds/"
query_directory = "/home/sareh/2020/sequences/ncbi/ncbi_pruned_genome/"
accession_file = "/home/sareh/2020/sequences/ncbi/test/test_accn.txt"
"""

# mpi
from mpi4py import MPI
the_world = MPI.COMM_WORLD
my_number = the_world.Get_rank()
total_number = the_world.Get_size()

complement_dict = {'A':'T', 'C':'G', 'G':'C', 'T':'A',
                    'W':'S', 'R':'Y', 'K':'M', 'Y':'R', 'S':'W', 'M':'K',
                    'B':'V', 'D':'H', 'H':'D', 'V':'B',
                    '*':'*', 'N':'N', '-':'-'}

def start_frame_shift(start,stop):

    if start-1 % 3 != 0:
        s = start - 2
        if s % 3 != 0:
            s = s - 1
        new_start = s
    else:
        new_start = start

    if stop + 1 %3 != 0:
        s = stop + 1
        if s + 1 %3 != 0:
            s = s + 1
        new_stop = s
    else:
        new_stop = stop

    return (int(new_start), int(new_stop))

def find_ovlp(ref_gene_directory):
    """
    run through all ref genes and their locations
    get a list of all locations and identifiy duplicates as ovlp regions
    more_ovlp is a list of indexes to remove to avoid frame shift of ovlp regions
    """
    all_loc = []
    for file in os.listdir(ref_gene_directory):
        ref = os.path.join(ref_gene_directory, file)
        ref_loc = parse_ref_header(ref)
        for i in ref_loc:
            a = i.split(":")
            start = int(a[0])
            stop = int(a[1])
        new_start = start
        new_stop = stop
        for n in range(int(new_start), int(new_stop)+1):
            all_loc.append(n)

    visited = set()
    dup = [x for x in all_loc if x in visited or (visited.add(x) or False)]
    unique_dup = set(dup)

    return(unique_dup)

def ref_ovlp_index(ref_loc, dup):
    """
    return list of index with no_ovlp regions, 0 index the cds
    """
    ref_all = []
    for i in ref_loc:
        a = i.split(":")
        start = a[0]
        stop = a[1]
        for n in range(int(start), int(stop)+1):
            ref_all.append(n)

    # Alternative
    no_ovlp=[]
    remove = len(ref_all) % 3
    end = len(ref_all) - remove + 1

    for index in range(3,end,3):
        a = index - 3
        b = index - 2
        c = index - 1

        if ref_all[a] in dup or ref_all[b] in dup or ref_all[c] in dup:
            continue
        #add = [a,b,c].sort()
        else:
            no_ovlp.extend([a,b,c])
    no_ovlp.sort()
    return(no_ovlp)

def ovlp(handle):
    all=[]
    #skip header [accn,prod1,loc1,dir1,prod2,loc2,dir2,seqlen1,seqlen2,overlap,shift]
    h = handle.readline().split(",")
    for line in handle:
        lst = (line.strip("\n")).split(",")
        all.append(lst)
    return(all)

def parse_ref_header(ref):
    h, original_rgene = get_ref(ref)
    header_lst = h.split(",")
    #>"NC_005300,NP_950235.1,glycoprotein precursor,1,92:5147"
    # ['NC_005300', 'NP_950235.1', 'glycoprotein precursor', '1', '92:5147']
    #>"NC_005219,NP_941987.1,envelope glycoprotein G2 (see comment),1,[94:1801](+)"
    # ['NC_005219', 'NP_941988.1', 'envelope glycoprotein G2 (see comment)', '1', '[1984:3418](+)']
    ref_loc = header_lst[-1].strip('"')
    lst = ref_loc.split(";") ##567:1147;1231:1440
    clean_lst = []
    for i in lst:
         a = re.search("[0-9]+\:[0-9]+",i) #to deal with ['[1984:3418](+)']
         b = a.group()
         clean_lst.append(b)

    return(clean_lst)

def reverse_and_complement(seq):
    rseq = seq[::-1]
    rcseq = ''
    for i in rseq:  # reverse order
        rcseq += complement_dict[i]
    return rcseq

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
    sequence = ""
    h = ""
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

def get_ref(ref):
    ref_handle = open(ref)
    for line in ref_handle:
        if line.startswith('>'):
            h = line.strip('>#\n')
        else:
            seq = line.strip('\n').upper()
    return(h, seq)

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
        result = query[rpos:(rpos+len(aligned))]
        if int(flag) & 16 != 0:
            # reverse-complement
            result = reverse_and_complement(result)
            #result = revcomp(result)
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


def cutter_minimap(ref, queries, outfile, out_ovlp, csvfile, no_ovlp_index):
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

    csvfile.write('header, p, original_refgene, len_nogap, len_noovlp, P_noovlp\n')

    # parse query genomes from input FASTA
    # queries = convert_fasta(fasta) #move to main

    # the reference cds seqeunce
    h, original_rgene = get_ref(ref)

    for qh, qs in queries:
        # minimap2 alignment: use minimap2 to extract homologous gene in query genome
        qgene_minimap = minimap2(query=qs, refseq=original_rgene)
        if len(qs) < 20:
            print("TOO SHORT qs " + str(len(qs)))
            print(qh)

        if qgene_minimap is None:
            # failed to align reference gene to query genome - no homology?
            # range = len(qs) * 0.1
            # trimmed_query = qs[int(start)-int(range):int(stop)+int(range)]
            trimmed_query = qs

            # if minimap fails skip minimap and do mafft pairwise alignment
            q1, r1 = mafft(query=trimmed_query, ref=original_rgene, trim=True)
            p1 = pdist(q1, r1)

            # check reverse-complement
            q2, r2 = mafft(query=reverse_and_complement(trimmed_query), ref=original_rgene, trim=True)
            p2 = pdist(q2, r2) #proportion of nucleotide differences, higher the more different

            if p1 < p2:
                qgene_mafft = q1
                rgene_mafft = r1
                p = p1
            else:
                qgene_mafft = q2
                rgene_mafft = r2
                p = p2

        else:
            # pairwise alignment
            qgene_mafft, rgene_mafft = mafft(query=qgene_minimap, ref=original_rgene, trim=False)
            p = pdist(qgene_mafft, rgene_mafft)
            sys.stdout.write("minimap: {} {} {}\n".format(len(rgene_mafft),len(qgene_mafft), qh))

        # to avoid out of index error of qgene_nogap += qgene_mafft[i]
        if len(qgene_mafft) != len(rgene_mafft):
            continue

        #if P > 0.5:  # P-distance threshold
        #    continue

        # remove gaps inserted into ref
        qgene_nogap = ""
        for i, rn in enumerate(rgene_mafft):
            if rn == '-': # skip insertion relative to reference
                continue
            qgene_nogap += qgene_mafft[i]

        # remove seq with >50% gap
        gap_per = qgene_nogap.count("-")/len(qgene_nogap)*100

        if gap_per < 50:
            outfile.write(">{}\n{}\n".format(qh, qgene_nogap))

            # remove ovlp
            qgene_no_ovlp = ""
            for i, rn in enumerate(qgene_nogap):
                if i in no_ovlp_index:
                    qgene_no_ovlp += rn
            out_ovlp.write(">{}\n{}\n".format(qh, qgene_no_ovlp))

            rgene_noovlp = ""
            for i, rn in enumerate(rgene_mafft):
                if i in no_ovlp_index:
                    rgene_noovlp += rn

            # Print length of each gene in the gapped outfile
            gap_qgene_mafft = qgene_mafft.count("-")
            gap_rgene_mafft = rgene_mafft.count("-")

            # calculate p-distance as a measure of similarity/alignment score
            p = pdist(qgene_mafft, rgene_mafft)
            p_noovlp = pdist(qgene_no_ovlp, rgene_noovlp)

            csvfile.write('{},{},{:1.3f},{},{},{},{:1.3f}\n'.format(qh, gap_per, 100*p,len(original_rgene), len(qgene_nogap),len(qgene_no_ovlp), 100*p_noovlp))

def main():
    count = 0
    accns = get_accn(accession_file) # all of the files to examine

    for accn in accns:
        print(accns)
        ref_gene_directory = os.path.join(ref_home_directory, accn)

        # Path to FASTA with genome sequences to process
        fasta = os.path.join(query_directory, accn)

        # if the ref_cds directory of the accn dosen't exsist skip
        if os.path.exists(ref_gene_directory) == False or os.path.exists(fasta) == False:
            continue

        # examine function find_ovlp (to print steps out)
        for file in os.listdir(ref_gene_directory):
            ref = os.path.join(ref_gene_directory, file)
            h, original_rgene = get_ref(ref)
            ref_loc = parse_ref_header(ref)

        dup = find_ovlp(ref_gene_directory)

        # FASTA with genome sequences to process
        fasta_handle = open(fasta)
        queries = convert_fasta(fasta_handle)

        cds_count = 0
        for file in os.listdir(ref_gene_directory):

            try:
                count += 1
                if count % total_number != my_number:
                    continue

                outfile = os.path.join(outdir, file)
                if os.path.isfile(outfile):
                    # output file exists, skip to next job
                    continue

                # Path to Reference CDS
                ref = os.path.join(ref_gene_directory, file)
                ref_handle = open(ref)

                # ref CDS coordinates (417:1546;1545:2132)
                ref_loc = parse_ref_header(ref)
                # 0 index of ref CDS without ovlp [0,1,2,7,8,9]
                no_ovlp_index = ref_ovlp_index(ref_loc, dup)

                # remove genes too short after cutting ovlp
                if len(no_ovlp_index) < 50:
                    continue
                cds_count += 1

                # Path to write FASTA output (cut gene sequences)
                outfile_handle = open(outfile, 'w+')  # append

                # Path to write FASTA output (NO OVLP)
                n = file+".noovlp"
                out_ovlp = os.path.join(outdir, n)
                out_ovlp_handle = open(out_ovlp, 'w+')  # append

                # Path to write CSV output (alignment scores)
                a = file+".csv"
                csvfile = os.path.join(outdir, a)
                csvfile_handle = open(csvfile, 'w+')

                # progress monitoring
                sys.stdout.write("[{} {}/{}] starting job {} {}\n".format(datetime.now().isoformat(), my_number, total_number, count, file))
                sys.stdout.flush()

                # parse query genomes from input FASTA
                # queries = convert_fasta(fasta_handle)
                cutter_minimap(ref, queries, outfile_handle, out_ovlp_handle, csvfile_handle, no_ovlp_index)

                ref_handle.close()
                outfile_handle.close()
                out_ovlp_handle.close()
                csvfile_handle.close()

            except:
                continue 


if __name__ == "__main__":
    main()
