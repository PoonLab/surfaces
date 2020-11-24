import os
import csv
from Bio import Entrez, SeqIO, SeqFeature
from csv import DictWriter
from time import sleep

Entrez.email = 'sbagher4@uwo.ca'
directory = '/Users/sarehchimeh/Data/Neighbours'

def parse_table(tbl):
    """
    Yield every row on the table parsed for downstream analysis 
    """
    handle = open(tbl)
    skip_line = next(handle)
    header = next(handle).strip('\n').split('\t')
    keywords = list(map(lambda x: x.strip('"'), header))
    print(keywords)

    row = None
    accno = ''

    for line in handle:
        delimiter = '\t' if '\t' in line else '     '
        values = line.strip('\n').split(delimiter)

        #Prepare next entry
        row = dict(zip(keywords, values))
        yield(row)


def retrieve_gid(accn):
    """
    Get genome ID for an accession number
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
    try:
        handle = Entrez.efetch(db='nucleotide', rettype='gb', retmode='text', id=gid)
        gb = SeqIO.parse(handle, format='genbank')
        return next(gb)
    except:
        print("Failed, trying again")
        return retrieve_record(gid)

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

def retrieve_annotations(record):
    an = record.annotations
    an_dict = {'Length': len(record.seq) , 'Topology':an['topology'],
                'Taxonomy':an['taxonomy'], 'Molecule': an['molecule_type'],
                'Proteins': len(record.features)
                }
    return an_dict

def main():
    virus_count = 0
    row_count = 0
    accn_count = 0
    refsecs = []
    for filename in os.listdir(directory):
        if filename.startswith('NC_') and len(filename) == 13:
            virus_count+=1
            file_path = '/Users/sarehchimeh/Data/Neighbours/{}'.format(filename)
            handle_one = open('/Users/sarehchimeh/Data/scraper/virus_info_{}.csv'.format(filename), 'w+')
            out_file_1 = DictWriter(handle_one, delimiter=',',
                                 extrasaction='ignore',
                                 fieldnames=[
                                     'Representative', 'Neighbor', 'Host', 'Selected lineage',
                                     'Taxonomy name', 'Length', 'Proteins',
                                     'Topology', 'Taxonomy', 'Molecule', 'Segment name"i'
                                 ])
            out_file_1.writeheader()

            handle_two = open('/Users/sarehchimeh/Data/scraper/orfs_location_{}.csv'.format(filename), 'w+')
            handle_two.write('accno,product,strand,coords,start_codon\n')

            for row in parse_table(file_path):
                accnos = [row['Neighbor']]
                row_count+=1
                print('row count is {}'.format(row_count))

                for accn in accnos:
                    accn_count+=1
                    print('accn count is {}'.format(accn_count))

                    if accn not in refsecs:
                        # only retrieves unique seqs
                        gid = retrieve_gid(accnos)

                        if gid is None:
                            print('Warning, failed to retrieve gid for {}'.format(accn))
                            continue

                        print(row['Taxonomy name'], gid, accn)
                        sleep(0.5)

                        record = retrieve_record(gid)
                        sleep(1)
                        annotation = retrieve_annotations(record)
                        sleep(0.5)
                        final_row = {**row, **annotation}
                        out_file_1.writerow(final_row)
                        handle_one.flush()
                        refsecs.append(accn)

                        for locus, product, strand, parts, start_codon, aaseq in retrieve_CDS(record):
                            parts_str = ';'.join('{}:{}'.format(p[0].position, p[1].position) for p in parts)
                            handle_two.write('{},"{}",{},{},{}\n'.format(accn, product[0], strand, parts_str, start_codon))

                        sleep(0.5)
    print('virus count is {}'.format(virus_count))

if __name__ == "__main__":
    main()
