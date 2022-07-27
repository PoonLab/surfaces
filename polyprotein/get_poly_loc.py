
from Bio import SeqIO
import os
import argparse
import csv

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--loc', type=str,
                       help='csv with location of polyprotiens')
    parser.add_argument('--csv', type=str,
                       help='csv identifying segmented and polyproteins')
    parser.add_argument('--dir', type=str,
                       help='directory of genbank record files')
    return parser.parse_args()

args = parse_args()

for file in os.listdir(args.dir):
    file_path = os.path.join(args.dir,file)
    record = SeqIO.read(file_path, "gb")

    line = []
    d = {}
    
    line.append(file)
    d["id"] = file

    for feature in record.features:
        # check if it's a segmented virus
        if feature.type == "source" and 'segment' in feature.qualifiers:
            line.append("seg")
            d["seg"] = "TRUE"
            #print(feature.qualifiers['segment'])
        else:
            line.append("nonseg")
            d["seg"] = "FALSE"

    for feature in record.features:
        # check is it's a polyprotein
        if feature.type == 'mat_peptide':
            d["poly"] = "TRUE"
            
            try:
                loc = feature.location
                prot_accn = feature.qualifiers['protein_id'][0]
                line = "{},{},{}\n".format(record.id, prot_accn, loc)
                loc_handle = open(args.loc, 'a')
                loc_handle.write(line)

            except:
                # skip viruses with missing record ID 
                print("issue: " + record.id)
        else: 
            d["poly"] = "FALSE"

    #write out the dictionary 
    l = ("{},{},{}\n".format(d["id"], d["seg"], d["poly"]))
    csv_handle = open(args.csv, "a")
    csv_handle.write(l)

print(d)
