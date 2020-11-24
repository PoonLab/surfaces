"""
retreieve protein coding sequences in fasta form
all into one output file
"""

from Bio import Entrez, SeqIO, SeqFeature
from csv import DictWriter
from time import sleep

Entrez.email = 'sbagher4@uwo.ca'

#Table with accession numbers
table = '/Users/sarehchimeh/Data/filtered_Neighbours.txt'

def pase_filtered_neighbour_table(table):
    """
    #:param table: filtered_Neighbours.txt
    #:return: list of all the neighbour accessions
    """
    lst=[]
    with open (table, 'r') as f:
        title = next(f)
        header = next(f)
        reader=csv.reader(f,delimiter='\t')
        for line in reader:
            lst.append(line[1])

    return(lst)

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


def retrieve_CDS(record):
    """
    Analyze features in Genbank record to extract (1) the number of coding
    regions (CDS)
    """
    cds = [feat for feat in record.features if feat.type=='CDS']
    for cd in cds:
        q = cd.qualifiers
        parts = []
        for part in cd.location.parts:
            parts.append((part.start, part.end))
        locus = q.get('locus_tag', '')
        product = q.get('product', [''])
        aaseq = q.get('translation', [''])
        yield locus, product, cd.strand, parts, q['codon_start'], aaseq


def main():

    # Protein file
    orffile = open('/Users/sarehchimeh/Data/lara_protein_scraper.fasta', 'w')
    #orffile.write('name,accno,product,strand,coords,aaseq\n')
    values = pase_filtered_neighbour_table(table)

    for accn in values:

        gid = retrieve_gid(accn)
        if gid is None:
            print('Warning, failed to retrieve gid for {}'.format(accn))
            continue

        print(accn)
        sleep(2)  # avoid spamming the server

        record = retrieve_record(gid)

        for locus, product, strand, parts, start, aaseq in retrieve_CDS(record):
            parts_str = ';'.join('{}:{}'.format(p[0].position, p[1].position) for p in parts)
            name = record.annotations['organism']
            product = product[0]
            seq = seq[0]
            orffile.write(f'>"{name},{accn},{product},{strand},{parts_str}"\n{aaseq}\n')

        sleep(2)


if __name__ == "__main__":
    main()