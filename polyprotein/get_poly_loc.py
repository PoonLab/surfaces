from Bio import SeqIO
import os
import argparse

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--out', type=str,
                       help='text files of genbank record')
    parser.add_argument('--dir', type=str,
                       help='directory')
    return parser.parse_args()

#dir = "/home/sareh/2020/polyprotein/gb_record/"
#outpath = "/home/sareh/2020/polyprotein/poly_loc2.csv"

args = parse_args()

for file in os.listdir(args.dir):
  file_path = os.path.join(args.dir,file)

  record = SeqIO.read(file_path, "gb")

  count = 0
  for feature in record.features:
     if feature.type == 'mat_peptide':
        count += 1
        try:
           loc = feature.location
           prot_accn = feature.qualifiers['protein_id'][0]
           #print(feature)
        except:
           print("issue")
           print(record.id)
           print(feature)
        if count > 1:
           outfile = open(args.out, 'a')
           line = "{},{},{}\n".format(record.id, prot_accn, loc)
           outfile.write(line)
