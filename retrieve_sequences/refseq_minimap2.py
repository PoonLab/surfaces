"""
retrieve full genome sequence of the refseq for each genome
each virus specie in it's own file "refseq_accn.fasta"
"""
import os
import csv
from Bio import Entrez, SeqIO, SeqFeature
from csv import DictWriter
from time import sleep

Entrez.email = 'sbagher4@uwo.ca'
ref_seq = '/Users/sarehchimeh/Data/Human_virus_accession.txt'

#from filter_neighbours.py
def parse_accesion(file):
    # Loop trough table to retrieve accession numbers as a list
    handle = open(file)
    accn = []
    for line in handle:
        temp = line.strip('\n')  # Remove new line
        accn.append(temp)
    return (accn)

#from full_genome_seq.py
def retrieve_gid(accn):
    """
    Find records according to genbank accession number
    """
    handle = Entrez.esearch(db='nucleotide', term=accn, retmax=1)
    response = Entrez.read(handle)
    if response['Count'] == '0':
        # retry query
        handle = Entrez.esearch(db='nucleotide', term=accn, retmax=1)
        response = Entrez.read(handle)
        if response['Count'] == '0':
            return None

    return response['IdList'][0]

#from full_genome_seq.py
def retrieve_record(gid):
    handle = Entrez.efetch(db='nucleotide', rettype='fasta', retmode='text', id=gid)
    gb = SeqIO.parse(handle, format='fasta')
    return next(gb)


def main():
    ref_seq_lst = parse_accesion(ref_seq)
    for accn in ref_seq_lst:
        with open('/Users/sarehchimeh/Data/refseq_genome/refseq_{}'.format(accn), 'w') as outfile:
            print(outfile)
            gid = retrieve_gid(accn)
            if gid is None:
                print('failed to retrieve gid for {}'.format(accn))
                continue
            sleep(1)
            record = retrieve_record(gid)
            SeqIO.write(record, outfile, 'fasta')
        sleep(2)


if __name__ == "__main__":
    main()
