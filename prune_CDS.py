#! /usr/bin/python3
#prune a fasta file from (Data/CDS nucleotide sequence)
# to keep only the sequence entries of the accn sequences in Data/pruned_accn
# input Data/CDS nucleotide sequence & Data/pruned_accn directories

import os
import re
from Bio import SeqIO

def get_accn(file):
    """
    Stream through txt file retrieve accns into a list
    """
    lst_accns = []
    with open(file, 'r') as f:
        for line in f:
            line = line.strip('\n')
            lst_accns.append(line)
            a = tuple(lst_accns)
        return(a)

# if "accn" in tuple:

def main():
    accn_directory = '/home/sareh/data/pruned_accns2'
    cds_directory = '/home/sareh/data/CDS_nucleotide_seq2'
    for filename1 in os.listdir(accn_directory):
        print('filename1 is '+ filename1)
        #filename1='pruned_accn_NC_007605'
        accn1=filename1[12:-3]
        print('accn1 ' + accn1)
        path1 = '{}/{}'.format(accn_directory, filename1)
        path2 = '{}/{}.fasta'.format(cds_directory, accn1)
        print(path2)
        # filename2='NC_001526.fasta'
        accessions= get_accn(path1)
        print(accessions)
        #all accessions in a tuple
        with open('/Users/sarehchimeh/Desktop/Pruned_CDS/Pruned_CDS_{}'.format(accn1),'w') as out_file:
            pruned_record_count = 0
            record_count = 0
            with open(path2,'r') as f2:
                for record in SeqIO.parse(path2, 'fasta'):
                    record_count += 1
                    id = []
                    a = record.id
                    x = re.findall('[A-Z]{1,3}[0-9]{5,6}.[0-9]', a)
                    # accns in a list
                    for i in x:
                        if i in accessions:
                            pruned_record_count +=1
                            SeqIO.write(record,out_file,'fasta')
        print('record count: '+ str(record_count))
        print('pruned record count: '+ str(pruned_record_count))
                                    # write the next line into a file (record sequence)

if __name__ == "__main__":
    main()
