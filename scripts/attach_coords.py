from Bio import SeqIO
import csv
import os
from glob import glob
from Bioplus.consensus import conseq
from Bioplus.pair_align import mafft
import tempfile
import sys


# parse annotation CSV
accns = {}
handle = open("/home/paula/surfaces_project/repos/surfaces/data/annotation_with_uniprot.csv", 'r')
reader = csv.DictReader(handle)
for row in reader:
    virus = row['virus']
    if virus not in accns:
        accns.update({virus: {}})
    accn = row['refseq'].split('.')[0]
    accns[virus].update({row['protein']: accn})
handle.close()


# retrieve sequences from Genbank files
refseqs = {}
for fn in glob("/home/surfaces/1_rawdata/refseq/*.gb"):
    record = SeqIO.read(fn, 'genbank')
    accn = os.path.basename(fn).split('.')[0]
    refseqs.update({accn: record.seq})


# process fasta files
types = ("*.fasta", "*.fa")
path5 = r"/home/surfaces/5_alt/"
path4 = r"/home/surfaces/4_filter/"

# retrieve step 5 files
files_5 = []
for t in types:
    files_5.extend(glob(f"{path5}{t}"))

files_4 = []
for t in types:
    files_4.extend(glob(f"{path4}{t}"))

virus_5 = set()
for fn in files_5:
    virus_5.add(os.path.basename(fn).split('_')[0])

virus_4 = set()
for fn in files_4:
    virus_4.add(os.path.basename(fn).split('_')[0])

missing = virus_4 - virus_5
print(f"these viruses were present in step4 folder but not in step5 folder: \n {missing}")

# print(f"there are {len(files_5)} files")

added = 0
for fn in files_4:
    if os.path.basename(fn).split('_')[0] in missing:
        # print(fn)
        added +=1
        files_5.append(fn)
print(f"there are {len(files_5)} files with an addition of {added}")

coord_dict = {}
for fn in files_5:
    # map to accession
    tokens = os.path.basename(fn).split('_')
    virus = tokens[0]
    if virus not in accns: # ask about potato versus potatoY virus are they the same?
        sys.stderr.write(f"ERROR: Could not find virus {virus} in annotations\n")
        sys.exit()
    
    protein = ' '.join(tokens[1:-1])
    if protein not in accns[virus]:
        sys.stderr.write(
            f"ERROR: Could not find {virus} protein '{protein}' in annotations\n")
        sys.exit()

    accn = accns[virus][protein]

    # print(f"processing {virus} {protein}")
    coord_dict.setdefault(virus, {})
    coord_dict[virus].setdefault(protein, {})


    # retrieve reference
    refseq = refseqs.get(accn, None)
    if refseq is None:
        sys.stderr.write(f"WARNING: failed to retrieve reference for {accn}\n")
        continue

    # generate consensus sequence
    aln = SeqIO.parse(fn, 'fasta')
    cseq = conseq(aln, thresh=0.5, is_nuc=True) # Bioplus
      
    # align to corresponding reference
    aquery, aref = mafft(cseq, refseq)
    gap = 0
    
    assert len(aquery) == len(aref) # ensure lengths match up
    # print("query and ref are the same length")

    ref_pos = 1
    start = None
    ranges = []
    for i, nuc in enumerate(aquery):
        if nuc != '-':
            if start is None:
                start = ref_pos
            ref_pos += 1
        else:
            if start is not None:
                stop = ref_pos - 1
                ranges.append((start, stop))
                start = None
            ref_pos += 1
    # handle final pos
    if start is not None:
        stop = ref_pos - 1
        ranges.append((start, stop))

    mapped_start = ranges[0][0]
    mapped_stop = ranges[-1][-1]

    coord_dict[virus][protein].update({"mapped_start": mapped_start, "mapped_stop": mapped_stop, "ranges": ranges})

# print(coord_dict)
handle = open("/home/paula/surfaces_project/repos/surfaces/data/annotation_with_uniprot.csv", 'r')
outfile = open("/home/paula/surfaces_project/repos/surfaces/data/annotation_with_uniprot_coords.csv", 'w')
reader = csv.DictReader(handle)
header = reader.fieldnames
header.extend(['mapped_start', 'mapped_stop', 'ranges'])

print(header)
writer = csv.DictWriter(outfile, fieldnames=header)
writer.writeheader()
for row in reader:
    virus = row['virus']
    protein = row['protein']
    coord_data = coord_dict[virus][protein]
    row.update({"mapped_start": coord_data["mapped_start"], "mapped_stop": coord_data["mapped_stop"], "ranges": coord_data["ranges"]})

    writer.writerow(row)
handle.close()
outfile.close()