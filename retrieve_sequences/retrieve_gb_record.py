import os
import re
import csv
from Bio import Entrez, SeqIO, SeqFeature
from csv import DictWriter
from time import sleep
from glob import glob
import argparse

Entrez.email = 'sbagher4@uwo.ca'
directory = '/Users/sarehchimeh/Data/Neighbours'

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--accns', type=str,
                       help='text files of accession numbers')
    return parser.parse_args()

def get_accns(file):
    """
    #:param table: filtered_Neighbours.txt
    #:return: list of all the neighbour accessions
    ex.['KF853231', 'DQ479956', 'JN704703', 'KF558382']
    """
    lst = []
    f = open(file, 'r')
    for line in f:
        name = line.split("\n")
        lst.append(name[0])
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
    handle = Entrez.efetch(db="nuccore", id=gid, rettype="gb", retmode="text")
    gb = SeqIO.parse(handle, format='genbank')
    return next(gb)

def main():
    args = parse_args()
    values = get_accns(args.accns)

    for accn in values:
       print(accn)

       # retrieve gb ID
       gid = retrieve_gid(accn)
       print(gid)
       if gid is None:
          print('Warning, failed to retrieve gid for {}'.format(accn))
          continue
          print(accn + ' ' + gid)
          sleep(2)

       # retrieving the gb record 
       handle  = Entrez.efetch(db="nuccore", id=gid, rettype="gb",retmode="text")
       text_gb = handle.read()

       out_path = os.path.join("./gb_record",accn)

       out_file = open(out_path ,"w")
       out_file.write(text_gb)
       out_file.close()

if __name__ == "__main__":
    main()
