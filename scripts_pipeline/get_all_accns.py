# Retrieve medatada, CDSs, and full genome from accession numbers
# Modified from: https://github.com/PoonLab/ProtClust/blob/main/scripts/protein_scraper.py#L100
from Bio import Entrez, SeqIO, SeqFeature
from time import sleep
import csv
import sys

import argparse


def get_args(parser):
    # Required arguments
    parser.add_argument(
        'table',
        help='Path to the file containing list of accession numbers'
    )
    parser.add_argument(
        'email',
        help='email name to be identified in NCBI'
    )
    # Optional
    parser.add_argument(
        '--outfile', default=None, 
        help='<Optional> Path to fasta file with amino acid sequences'
    )
    parser.add_argument(
        '--poly', action = 'store_true', default = False, 
        help='<Optional> true if virus encodes a polyprotein'
    )

    return parser.parse_args()

def parse_table(tbl):
    """
    Loop trough table file to retrieve accession numbers as a list
    :param tbl: list with accession numbers
    """
    handle = open(tbl)

    values = []
    for line in handle:
        temp = line.strip('\n')  # Remove new line
        values.append(temp)

    return(values)

def retrieve_gid(accn):
    """
    Find records according to genbank accession number
    :param accn: NCBI accession number
    :returns: GenInfo Identifier (GI) associated with the accession number
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
    """
    Using Entrez, fetch SeqRecord object from NCBI using GI number
    :param gid: GenInfo Identifier (GI) assigned to a record proccesed by NCBI 
    :returns: SeqRecord object 
    """
    handle = Entrez.efetch(db='nucleotide', rettype='gb', retmode='text', id=gid)
    # handle = Entrez.efetch(db='nuccore', id=gid, rettype='gb', retmode='text')
    gb = SeqIO.parse(handle, format='genbank')
    return next(gb)

def retrieve_CDS(record, poly):
    """
    Analyze features in Genbank record to extract (1) the number of coding regions (CDS)
    :param record: SeqRecord object (used in BioPython to hold a sequence and sequence information)
    """

    if poly:
        cds = [feat for feat in record.features if feat.type=="mat_peptide"] 
    else:
        cds = [feat for feat in record.features if feat.type=='CDS']
    for cd in cds:
        q = cd.qualifiers
        parts = []
        for part in cd.location.parts:
            parts.append((part.start, part.end))
        locus = q.get('locus_tag', '')
        product = q.get('product', [''])
        cd_seq = cd.location.extract(record).seq
        if cd_seq[0:3] == "ATG" and not poly:
            aaseq = cd_seq.translate(cds=True)
            yield locus, product, cd.strand, parts, cd_seq, aaseq
        else:
            aaseq = cd_seq.translate(cds=False)
            stop_count = aaseq.count("*")
            if stop_count > 3:  # More than three stop codons in the amino acid seq
                print(f"\nSkipping record {record.name}")
                print(f"CDS with stop codons:")
                print(aaseq, "\n")
                break
            else:
                yield locus, product, cd.strand, parts, cd_seq, aaseq

def get_metadata(accn):
    handle = Entrez.efetch(db='nuccore', id=accn, rettype='gb', retmode='text')
    records = SeqIO.parse(handle, 'gb')
    
    for rec in records:
        taxon = ':'.join(rec.annotations["taxonomy"])
        source = list(filter(lambda f: f.type=='source', rec.features))
        quals = source[0].qualifiers
        iso = quals.get('isolate', [''])[0],
        host = quals.get('host', [''])[0],
        loc = quals.get('country', [''])[0],
        coldate = quals.get('collection_date', [''])[0]
        organism = quals.get('organism', [''])[0]
        yield rec.id, taxon, iso, host, loc, coldate, rec.seq, organism


def main():

    parser = argparse.ArgumentParser(
        description='From list of accession numbers, retrieve amino acid sequences as fasta file'
    )
    
    args = get_args(parser)
    accn_table = args.table
    Entrez.email = args.email
    filename = args.outfile if args.outfile else f'{accn_table}'
    poly = args.poly
    
    cds_file = open(f'{filename}_CDSs.fasta', 'w')
    aa_file = open(f'{filename}_aa.fasta', 'w')
    md_file = open(f'{filename}_md.csv', 'w')
    md_writer = csv.writer(md_file)
    md_writer.writerow(['accn', 'taxonomy', 'isolate', 'host', 'country', 'coldate'])
    values = parse_table(accn_table)
    # full_gen_file = open(f'{filename}_full_gen.fasta', 'w')

    for accn in values:

        gid = retrieve_gid(accn)
        if gid is None:
            print('Warning, failed to retrieve gid for {}'.format(accn))
            continue

        sleep(1)  # avoid spamming the server
        record = retrieve_record(gid)
        
        # Get medatada and store full genomes
        for id, taxon, iso, host, loc, coldate, full_seq, organism in get_metadata(accn):
            md_writer.writerow([id, taxon, iso, host, loc, coldate])
            # full_gen_file.write(f'>{id}-{organism}-{host}-{coldate}\n{full_seq}\n')

        # get sequences
        for locus, product, strand, parts, cd_seq, aa_seq in retrieve_CDS(record, poly):
            location = '.'.join('{}_{}'.format(p[0], p[1]) for p in parts)
            organism_name = record.annotations['organism'].replace(" ", "_")
            product = product[0].replace(" ", "_")
            cds_file.write(f'>{accn}-{organism_name}-{product}-{strand}-{location}\n{cd_seq}\n')
            aa_file.write(f'>{accn}-{organism_name}-{product}-{strand}-{location}\n{aa_seq}\n')

        print(accn)
        sleep(1)

    print(f'\n>>> DONE <<<\n')

if __name__ == "__main__":
    main()
