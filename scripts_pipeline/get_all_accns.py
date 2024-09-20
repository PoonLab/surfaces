description = """
Retrieve medatada, coding sequences (CDSs), and full genome 
sequences from the NCBI Genbank database, given a list of accession 
numbers
"""

from Bio import Entrez, SeqIO
from time import sleep
import csv
import sys
import argparse
import re


def read_accns(handle, regex="([A-Z]+_?[0-9]+)(\.[0-9]+)?"):
    """ Read and validate accession numbers from file """
    pat = re.compile(regex)
    accns = []
    for line in handle:
        matches = pat.findall(line.strip())
        if len(matches) == 0:
            if line != '\n':  # ignore blank line
                sys.stderr.write(f"Invalid accession detected: {line}\n")
            continue
        accns.append(''.join(matches[0]))
    return accns


def get_metadata(record):
    """
    Retrieve source metadata associated with Genbank SeqRecord.
    :return:  dict, values associated with fields
    """    
    fields = ('isolate', 'host', 'country', 'collection_date', 'organism')
    source = record.features[0]
    values = [source.qualifiers.get(f, [''])[0] for f in fields]
    
    result = dict(zip(fields, values))
    result['accn'] = record.id
    result.update({'taxonomy': ':'.join(record.annotations["taxonomy"])})

    return result


def retrieve_CDS(record):
    """
    Analyze features in Genbank record to extract (1) the number of coding regions (CDS)
    :param record: SeqRecord object (used in BioPython to hold a sequence and sequence information)
    """
    cds = [feat for feat in record.features if feat.type=='CDS']
    for cd in cds:
        q = cd.qualifiers
        parts = []
        
        # return cds annotation
        for part in cd.location.parts:
            parts.append((part.start, part.end))
        product = q.get('product', [''])
        cd_seq = cd.location.extract(record).seq
        
        # translate aa sequence
        aaseq = cd_seq.translate(cds=False)
        stop_count = aaseq.count("*")
        if stop_count > 3:  # More than three stop codons in the amino acid seq
            print(f"\nSkipping record {record.name}")
            print(f"CDS with {stop_count} stop codons:")
            print(aaseq, "\n")
            break

        yield product, cd.strand, parts, cd_seq, aaseq


if __name__ == "__main__":
    # command line interface
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('infile', type=argparse.FileType('r'),
                        help='Path to file containing accession numbers')
    parser.add_argument('email', help='email for Entrez API transactions')
    parser.add_argument('--prefix', type=str, 
                        help='Path and prefix for output files')
    parser.add_argument("--batch", type=int, default=50,
                        help="Batch size for retrieving records")
    parser.add_argument("--delay", type=float, default=1.0, 
                        help="Number of seconds to pause between queries (default 1).")
    parser.add_argument('--poly', action='store_true',
                        help='Set if virus genome encodes a polyprotein')
    parser.add_argument('--original', action='store_true', 
                        help="Keep original sequence headers")
    args = parser.parse_args()
    
    Entrez.email = args.email

    # read accessions from file
    accns = read_accns(args.infile)
    sys.stderr.write(f"Extracted {len(accns)} accession numbers from input file:\n")
    sys.stderr.write(f"  {', '.join(accns[:3])}, ..., {', '.join(accns[-3:])}\n")
    
    # prepare output file paths
    if args.prefix is None:
        args.prefix = args.infile.name
    cds_file = open(f'{args.prefix}_CDSs.fasta', 'w')
    aa_file = open(f'{args.prefix}_aa.fasta', 'w')
    if args.poly:
        cds_file = open(f'{args.prefix}_CDSs_polyprot.fasta', 'w')
        aa_file = None

    fields = ['accn', 'taxonomy', 'isolate', 'host', 'country', 'collection_date']
    md_file = open(f'{args.prefix}_md.csv', 'w')
    md_writer = csv.DictWriter(md_file, fieldnames=fields)
    md_writer.writeheader()

    # retrieve records in batches
    for i in range(0, len(accns), args.batch):
        query = ','.join(accns[i:(i+args.batch)])
        response = Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", id=query)
        for record in SeqIO.parse(response, format="genbank"):
            mdata = get_metadata(record)
            org = mdata.pop('organism').replace(' ', '_')
            md_writer.writerow(mdata)
            
            # get sequences
            for prod, strand, parts, cd_seq, aa_seq in retrieve_CDS(record):
                loc = '.'.join('{}_{}'.format(p[0], p[1]) for p in parts)
                prod = prod[0].replace(" ", "_").replace("-", "_")
                if args.original:
                    cds_file.write(f">{record.description}\n{cd_seq}\n")
                    if not args.poly:
                        aa_file.write(f">{record.description}\n{aa_seq}\n")
                else:
                    cds_file.write(f'>{record.id}-{org}-{prod}-{strand}-{loc}\n{cd_seq}\n')
                    if not args.poly:
                        aa_file.write(f'>{record.id}-{org}-{prod}-{strand}-{loc}\n{aa_seq}\n')            

        sleep(args.delay)

    # close files
    cds_file.close()
    if not args.poly:
        aa_file.close()
    print('\n>>> DONE <<<\n')
