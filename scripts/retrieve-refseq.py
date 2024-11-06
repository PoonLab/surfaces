"""
Utility script for downloading RefSeq files
"""

from Bio import Entrez, SeqIO
import csv
import time


Entrez.email = 'apoon42@uwo.ca'
handle = open("data/annotation.csv")
reader = csv.DictReader(handle)

accns = set()
for row in reader:
    accns.add(row['refseq'])

for accn in accns:
    print(accn)
    response = Entrez.efetch(
        db='nucleotide', rettype='gb', retmode='text', id=accn)
    fn = f"data/refseq/{accn.split('.')[0]}.gb"
    with open(fn, 'w') as outfile:
        outfile.write(response.read())
    time.sleep(1)
