"""
retrieves coding sequences for the Accession
input directory: files with accessions
output directory: a directory is made for each accession with each CDS in a fasta file
"""
import os
import re
import csv
from Bio import Entrez, SeqIO, SeqFeature
from csv import DictWriter
from time import sleep
from glob import glob
import argparse

Entrez.email = 'sbagher4@uwo.ca'

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--accndir', default='/home/sareh/data/Human_virus_accession.txt', type=str,
                       help='file of accession numbers')
    parser.add_argument('--outdir', default='/home/sareh/data/nuc_refseq_CDS/', type=str,
                        help='Directory to write outputs.')
    return parser.parse_args()

def pase_filtered_neighbour_table(table):
    """
    #:param table: filtered_Neighbours.txt
    #:return: list of all the neighbour accessions
    ex.['KF853231', 'DQ479956', 'JN704703', 'KF558382']
    """
    lst = []
    with open(table, 'r') as f:
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
    handle = Entrez.efetch(db='nucleotide', rettype='gb', retmode='text', id=gid)
    gb = SeqIO.parse(handle, format='genbank')
    return next(gb)

def main():
    args = parse_args()
    # get all accns for a virus in a list
    values = pase_filtered_neighbour_table(args.accndir)
    for accn in values:
        path = ("{}/{}".format(args.outdir,accn))
        os.mkdir(path)
        print(accn)  # KF853231
        gid = retrieve_gid(accn)
        record = retrieve_record(gid)
       #retrieve each CDS
        cds = [feat for feat in record.features if feat.type == "CDS"]  # [SeqFeature(FeatureLocation(ExactPosition(154), ExactPosition(544), strand=1), type='CDS'), SeqFeature(FeatureLocation(ExactPosition(570), ExactPosition(897), strand=-1), type='CDS'), SeqFeature(FeatureLocation(ExactPosition(1115), ExactPosition(1832), strand=1), type='CDS')
        for cd in cds:
            q = cd.qualifiers  # OrderedDict([('gene', ['ORF 0']), ('codon_start', ['1']), ('product', ['ORF 0']), ('protein_id', ['AHJ09141.1']), ('translation', ['MATVHYSRRPGTPPVTLTSSPGMDDVATPIPYLPTYAEAVADAPPPYRSRESLVFSPPLFPHVENGTTQQSYDCLDCAYDGIHRLQLAFLRIRKCCVPAFLILFGILTLTAVVVAIVAVFPEEPPNSTT'])])
            proteinID = q.get('protein_id', [''])  # ['AHJ09212.1']
            protein_accn = proteinID[0]  # AHJ09141.1
            product = q.get('product', [''])[0]  # ORF 0
            locus = q.get('locus_tag', '')  # empty
            strand = cd.strand  # 1 or -1
            start = q['codon_start']  # ['1']
            seq = record.seq  # full nucleotide seqeunce
            name = record.annotations['organism']
            gene_sequence = cd.extract(record.seq)
            print(protein_accn)
            print(len(gene_sequence))
            parts = []  # [(ExactPosition(154), ExactPosition(544))]
            for part in cd.location.parts:  # [154:544](+)
                parts.append((part.start, part.end))
            parts_str = ';'.join('{}:{}'.format(p[0].position, p[1].position) for p in parts)

            with open('{}/{}_{}'.format(path, accn, protein_accn), 'w') as outfile:
                outfile.write(f'>"{accn},{protein_accn},{product},{strand},{parts_str}"\n{gene_sequence}\n')

if __name__ == "__main__":
    main()
