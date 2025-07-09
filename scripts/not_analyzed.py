import os
import csv
from glob import glob
from collections import defaultdict

types = ("*step3.fasta", "*step3.cds.fa", "measles*step2.fasta") # measles only step 2 fasta present
filepath = r"/home/surfaces/3_codaln/"
all_dict = defaultdict(list)
files_3 = []
for t in types:
    files_3.extend(glob(f"{filepath}{t}"))
    
for f in files_3:
    parts = os.path.basename(f).split('_')
    virus = parts[0]
    protein = ' '.join(parts[1:-1])
    all_dict[virus].append(protein)
print(all_dict)

fubar_filepath = "/home/surfaces/6_fubar"
fubar_files = glob(f"{fubar_filepath}/*.fubar.json")
fubar_dict = defaultdict(list)
for f in fubar_files:
    parts = os.path.basename(f).split('_')
    virus = parts[0]
    protein = ' '.join(parts[1:-1])
    fubar_dict[virus].append(protein)
print(fubar_dict)

dropped = defaultdict(list)
for virus, proteins in all_dict.items():
    for prot in proteins:
        if prot not in fubar_dict[virus]:
            print(type(prot))
            dropped[virus].append(prot)

print(dropped)

with open("not_analyzed.txt", 'w') as outfile:
    writer = csv.writer(outfile)
    writer.writerow(["virus", "protein"])
    for virus, proteins in dropped.items():
        for p in proteins:
            writer.writerow([virus, p])
