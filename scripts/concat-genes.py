import sys
import argparse
import re
from Bio import SeqIO


description = """
Influenza A/B virus encodes multiple open reading frames on some of its
segments.  For example, IAV segment 7 encodes M1 and M2, where M2 consists
of a small portion of the 5' end of M1, 
Concatenate alternative open reading frames and exclude overlaps.
For the current pipeline, sequence headers are delimited by underscores:
  {accn}-{organism}-{product}-{strand}-{location}
where location consists of closed intervals separated by '.' for spliced genes
"""

# regex for Genbank accession numbers
accpat = re.compile("([A-Z]{1,3}_?[0-9]{5,7})(\.[0-9]+)?")


def concat(records):
   # extract feature locations from sequence names
    coords = []
    for h, _ in records:
        accn, organism, product, strand, location = h.split('-')
        for loc in location.split('.'):
        
        coords.append([ for left, right in location.split('.')])
        matches = pat.findall(h)
        if len(matches) == 0:
            print(block)
            sys.exit()
        coords.append([int(i) for i in nums.findall(matches[0])])
    
    # find main (unspliced) sequence first (M1)
    finds = list(filter(lambda x: len(x)==2, coords))
    if len(finds) == 0:
        print("Error: failed to locate M1 sequence")
        print(records)
        sys.exit()
    
    m1 = finds[0]
    which_m1 = coords.index(m1)
    header, seq = records[which_m1]

    # append other sequence
    for i, co in enumerate(coords):
        if i == which_m1:
            continue
        nuclen = co[-1] - (m1[-1]+1)
        aalen = nuclen//3
        aaseq = records[i][1]
        seq += aaseq[-aalen:]
        break
    
    return header, seq



parser = argparse.ArgumentParser(
    description=description
)
parser.add_argument("infile", type=argparse.FileType('r'),
                    help="Path to file with CDS records")
parser.add_argument("-o", "--outfile", type=argparse.FileType('w'), 
                    default=sys.stdout,
                    help="Path to write output FASTA (default stdout)")
args = parser.parse_args()


# collect sequences in FASTA by accession
by_accn = {}
for record in SeqIO.parse(args.infile, 'fasta'):
    matches = accpat.findall(record.description)
    if len(matches) == 0:
        sys.stderr.write(f"ERROR: Failed to extract accession number for {record.description}\n")
        sys.exit()
    accn = ''.join(matches[0])
    if accn not in by_accn:
        by_accn.update({accn: []})
    by_accn[accn].append(record)

# concatenate sequences
for accn, seqs in by_accn.items():


sys.exit()

blocks = args.infile.read().split("\n\n")
for block in blocks:
    try:
        records = convert_fasta(block.split('\n'))
    except UnboundLocalError:
        if block == '':
            continue  # handle extra lines or end of file
        raise
    
    # trivial case, only one feature (ORF) in sequence
    if len(records) == 1:
        args.outfile.write(block+'\n')
        continue

    header, seq = concat(records)
    args.outfile.write(f">{header}\n{seq}\n")
