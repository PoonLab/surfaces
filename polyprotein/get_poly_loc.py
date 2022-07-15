from Bio import SeqIO
import os

dir = "/home/sareh/2020/polyprotein/gb_record/"
outpath = "/home/sareh/2020/polyprotein/poly_loc.csv"

for file in os.listdir(dir):
  file_path = os.path.join(dir,file)

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
           print(feature)
        if count > 1:
           outfile = open(outpath, 'a')
           line = "{},{},{}\n".format(record.id, prot_accn, loc)
           outfile.write(line)
