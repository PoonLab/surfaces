"""
retrieves full nucleotide sequences for the Accession
input file: all files from the Data/Nucleotide directory
output file: files in the Data/sequences Directory named nseq_NC_000001.txt
"""
import os
import csv
from Bio import Entrez, SeqIO, SeqFeature
from csv import DictWriter
from time import sleep

Entrez.email = 'sbagher4@uwo.ca'
directory = '/Users/sarehchimeh/Data/Neighbours'

def pase_filtered_neighbour_table(table):
    """
    #:param table: filtered_Neighbours.txt
    #:return: list of all the neighbour accessions
    ex.['KF853231', 'DQ479956', 'JN704703', 'KF558382']
    """
    lst = []
    with open(table, 'r') as f:
        title = next(f)
        header = next(f)
        reader = csv.reader(f, delimiter='\t')
        for line in reader:
            lst.append(line[1])

    return (lst)

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

def retrieve_record(gid):
    handle = Entrez.efetch(db='nucleotide', rettype='gb', retmode='text', id=gid)
    gb = SeqIO.parse(handle, format='genbank')
    return next(gb)

def main():
    for filename in os.listdir(directory):
        outfile = open('/Users/sarehchimeh/Data/sequences/nseq_{}'.format(filename), 'w')
        file_path = '/Users/sarehchimeh/Data/Neighbours/{}'.format(filename)
        values = pase_filtered_neighbour_table(file_path)

        for accn in values:
            gid = retrieve_gid(accn)
            if gid is None:
                print('Warning, failed to retrieve gid for {}'.format(accn))
                continue
                print(accn + ' ' + gid)
            sleep(2)  # avoid spamming the server
            record = retrieve_record(gid)
            name = record.annotations['organism']
            seq = record.seq
            outfile.write('>{},{},\n{}\n'.format(name, accn, seq)
        sleep(2)

if __name__ == "__main__":
    main()