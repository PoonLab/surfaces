#combine testminimap and cutter_mpi.py

from glob import glob
import tempfile
import subprocess
import re

#mpi
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

def convert_fasta (file):
    """
    takes in file path 
    """
    result = []
    sequence = ''
    with open(file,'r') as handle:
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


def cutter_minimap(ref, fasta, outfile, csvfile):
    """
    :param ref:  open stream in read mode to FASTA file containing reference
                 gene sequence
    :param fasta:  open stream in read mode to FASTA file containing query
                   genome sequences
    :param outfile:  open stream in write mode to output trimmed sequences
                     in FASTA format
    :param csvfile:  open stream in write mode to output alignment scores
    """
  
    scorefile.write('header,align.score\n')
    
    queries = convert_fasta(fasta)
    
    rgene = convert_fasta(ref)[0][1] #the reference cds seqeunce
    
    for qh, qs in queries:
        qgene = minimap2(query=qs, refseq=rgene) #each alignment 
        if qgene != None:
            out_file.write('>{}\n{}\n'.format(qh, qgene))
            try:
                ndiff = 0
                #print(rgene) #print(len(rgene)) #print(qgene) #print(len(qgene))
                p = pdist(qgene, rgene)
                scorefile.write('{},{:1.3f}\n'.format(qh, 100*p))
                #print('{},{:1.3f}'.format(qh, 100*p))
                worked_count += 1 
            except Exception as e:
                #print("i {}, nt1 {},nt2 {},qh {} \n {} \n {}".format(i,nt1,nt2,qh,qgene,rgene))
                sys.stdout.write("{},{}".format(qh,(len(rgene)-len(qgene))))
                #print("qgene length:{}".format(len(qgene))) #print("rgene length:{}".format(len(rgene)))
                #print("ERROR : "+str(e))
                error_count += 1
        else:
            sys.stdout.write("{}\n".format(qh))
            #print('{},{:1.3f}'.format(qh, 100*p))
            #print("{},{},{}".format(accn,worked_count,error_count))


def main():
    count = 0
    
    #Directories 
    outdir = '/home/sareh/surfaces/find_cds/minimap' #where all the files would be written into
    ref_home_directory = '/home/sareh/surfaces/find_cds/corrected_ref_cds'
    query_directory = '/home/sareh/data/pruned_genome'
    accession_file = '/home/sareh/surfaces/find_cds/virus_pruned_genome.txt'

    accns = get_accn(accession_file) #all of the files to examine

    for accn in accns:
        ref_gene_directory = os.path.join(ref_home_directory, accn)
        
        # Path to FASTA with genome sequences to process
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

            #os.path.isfile(path) | os.path.exists(path) | path.exists
            if os.path.isfile(outfile):
                # output file exists, skip to next job
                continue

            # progress monitoring
            sys.stdout.write("[{} {}/{}] starting job {} {}\n".format(
                datetime.now().isoformat(), my_number, total_number, count, file))
            sys.stdout.flush()

            # run analysis
            out_file = open(outfile, 'w+')  # append
            scorefile = open(csvfile, 'w+')
            ref_handle = open(ref)
            fasta_handle = open(fasta)
            cutter_minimap(ref, fasta, outfile, csvfile)
            #if os.path.isfile(outfile) and count % total_number == my_number:
            ref_handle.close()
            fasta_handle.close()
            out_file.close()
            scorefile.close()

def test():

    non_count = 0
    ran_count = 0
    all_count = 0
    error_count = 0 

    #Directories 
    ref_home_directory = '/home/sareh/surfaces/find_cds/corrected_ref_cds'
    query_directory = '/home/sareh/data/pruned_genome'
    accession_file = '/home/sareh/surfaces/find_cds/virus_pruned_genome.txt' 

    accns = get_accn(accession_file) #all of the files to examine
    score = '/home/sareh/surfaces/find_cds/cutter_mpi_minimap_score.txt'
    error = '/home/sareh/surfaces/find_cds/cutter_mpi_minimap_error.txt'
    non_query = '/home/sareh/surfaces/find_cds/cutter_mpi_minimap_none.txt'
    score_handle = open(score,'w')
    error_handle = open(error,'w')
    non_handle = open(non_query,'w')
    for accn in accns:
        ref_gene_files = glob('/home/sareh/surfaces/find_cds/corrected_ref_cds/{}/{}_*'.format(accn,accn)) 
        query_file_path = ('/home/sareh/data/pruned_genome/Pruned_nuc_{}'.format(accn))
        #query_file_handle = open(query_file_path)
        queries = convert_fasta(query_file_path)
        
        #Each gene files
        for rgf in ref_gene_files:
            #rfg_handle = open(rfg)
            rgene = convert_fasta(rgf)[0][1] #the reference cds seqeunce
            #print(rgene)
            #Each query header and sequence
            for qh, qs in queries:
                #print(qh)
                qgene = minimap2(query=qs, refseq=rgene) #each alignment 
                all_count += 1
                if qgene != None:
                    ran_count += 1
                    try:
                        ndiff = 0
                        #print(rgene) #print(len(rgene)) 
                        #print(qgene) #print(len(qgene))
                        p = pdist(qgene, rgene)
                        score_handle.write('{},{:1.3f}\n'.format(qh, 100*p))
                        #print('{},{:1.3f}'.format(qh, 100*p))
                    except Exception as e:
                        #print("i {}, nt1 {},nt2 {},qh {} \n {} \n {}".format(i,nt1,nt2,qh,qgene,rgene))
                        print("{},{}".format(qh,(len(rgene)-len(qgene))))
                        print("qgene length:{}".format(len(qgene)))
                        print("rgene length:{}".format(len(rgene)))
                        print("ERROR : "+str(e))
                        error_count += 1
                else:
                   non_handle.write("{}\n".format(qh))
                   non_count += 1 
                #print('{},{:1.3f}'.format(qh, 100*p))
            #print("{},{},{}".format(accn,worked_count,error_count))
    print("non count is {}".format(non_count))
    print("all count is {}".format(all_count))
    print("ran count is {}".format(ran_count))
    print("error count is {}".format(error_count))
if __name__ == "__main__":
    test()
