#! /usr/bin/python3
# prune a fasta file from (Data/CDS nucleotide sequence)
# to keep only the sequence entries of the accn sequences in Data/pruned_accn
# input Data/CDS nucleotide sequence & Data/pruned_accn directories

import os
import re
from Bio import SeqIO
from glob import glob
import argparse


def get_accn(file):
    """
    Stream through txt file retrieve accns into a list
    :param file:  str, path to file
    """
    lst_accns = []
    with open(file, 'r') as f:
        for line in f:
            line = line.strip('\n')
            lst_accns.append(line)
    a = tuple(lst_accns)
    return(a)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('accndir', default='/home/sareh/data/pruned_accns', type=str,
                       help='Directory containing text files of accession numbers from pruned trees.')
    parser.add_argument('cdsdir', default='/home/sareh/data/CDS_nucleotide_seq', type=str,
                        help='Directory containing FASTA files of protein coding sequences')
    parser.add_argument('outdir', default='/home/sareh/data/Pruned_CDS/', type=str,
                        help='Directory to write outputs.')
    return parser.parse_args()


def main():
    args = parse_args()
    accn_files = glob(os.path.join(args.accndir, 'pruned_accn_*'))
    
    # iterate through all files listing accessions grouped by virus species
    for path1 in accn_files:
        filename1 = os.path.basename(path1)
        print('filename1 is '+ filename1)
        #filename1='pruned_accn_NC_007605'
        
        # retrieve reference accession from file name
        accn1 = filename1.replace('pruned_accn', '')  # [12:]
        print('accn1 ' + accn1)
        
        #path1 = '{}/{}'.format(args.accndir, filename1)
        # use reference accn to retrieve corresponding CDS FASTA
        cds_path = os.path.join(args.cdsdir, accn1+'.txt')
        #path2 = '{}/{}.txt'.format(args.cdsdir, accn1)
        print(cds_path)
        
        # filename2='NC_001526.fasta'
        accessions = get_accn(path1)  # returns tuple of accessions in file
        print(accessions)
        
        # all accessions in a tuple
        accn_regex = re.compile('[A-Z]{1,3}[0-9]{5,6}.[0-9]')
        with open('{}/Pruned_CDS_{}'.format(args.outdir, accn1), 'w') as out_file:
            pruned_record_count = 0
            record_count = 0
            
            # construct dict from CDS FASTA
            cds_dict = {}
            with open(cds_path, 'r') as f2:
                for record in SeqIO.parse(cds_path, 'fasta'):
                    this_accn = accn_regex.findall(record.id)[0]
                    cds_dict.update({this_accn, str(record.seq)})
                
            

                for record in SeqIO.parse(path2, 'fasta'):
                    record_count += 1
                    id = []
                    this_id = record.id
                    this_accn = re.findall(, this_id)
                    # accns in a list
                    for i in x:
                        if i in accessions:
                            pruned_record_count +=1
                            SeqIO.write(record,out_file,'fasta')
        print('record count: '+ str(record_count))
        print('pruned record count: '+ str(pruned_record_count))
                                    # write the next line into a file (record sequenc$

if __name__ == "__main__":
    main()



