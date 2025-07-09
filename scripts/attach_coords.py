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
accns_set = set()
handle = open("/home/paula/surfaces_project/repos/surfaces/data/annotation_with_uniprot.csv", 'r')
reader = csv.DictReader(handle)
for row in reader:
    virus = row['virus']
    # if virus == "enterov1": # manully updated
    #     virus = "enteroA71"
    if virus not in accns:
        accns.update({virus: {}})
    accn = row['refseq'].split('.')[0]
    accns[virus].update({row['protein']: accn})
    accns_set.add(accn)
handle.close()


# retrieve sequences from Genbank files
refseqs = {}
for fn in glob("/home/surfaces/1_rawdata/refseq/*.gb"):
    record = SeqIO.read(fn, 'genbank')
    accn = os.path.basename(fn).split('.')[0]
    refseqs.update({accn: record.seq})


# # process fasta files
types = ("*step3.fasta", "*step3.cds.fa", "measles*step2.fasta") # measles only step 2 fasta present
filepath = r"/home/surfaces/3_codaln/"
files_3 = []
for t in types:
    files_3.extend(glob(f"{filepath}{t}"))

## deprecated - only using step 3 fasta 
# path5 = r"/home/surfaces/5_alt/"
# path4 = r"/home/surfaces/4_filter/"

# # retrieve step 4 & 5 files
# files_5 = []
# for t in types:
#     files_5.extend(glob(f"{path5}{t}"))

# files_4 = []
# for t in types:
#     files_4.extend(glob(f"{path4}{t}"))

# virus_5 = set()
# for fn in files_5:
#     virus_5.add(os.path.basename(fn).split('_')[0])

# virus_4 = set()
# for fn in files_4:
#     virus_4.add(os.path.basename(fn).split('_')[0])

# missing = virus_4 - virus_5
# print(f"these viruses were present in step4 folder but not in step5 folder: \n {missing}")

# # print(f"there are {len(files_5)} files")

# added = 0
# for fn in files_4:
#     if os.path.basename(fn).split('_')[0] in missing:
#         # print(fn)
#         added +=1
#         files_5.append(fn)
# print(f"there are {len(files_5)} files with an addition of {added}")


coord_dict = {}
aligns = {}
for fn in files_3:
    # map to accession
    tokens = os.path.basename(fn).split('_')
    virus = tokens[0]
    if virus not in accns: # ask about potato versus potatoY virus are they the same?
        sys.stderr.write(f"ERROR: Could not find virus {virus} in annotations\n")
        sys.exit()
    
    protein = ' '.join(tokens[1:-1])
    if protein not in accns[virus]:
        if ".fasta" in protein or ".fa" in protein: # handling misnamed files 
            protein = protein.split('.')[0]
        else:
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

    ref_seqid = f"ref_{virus}_full"
    query_seqid = f"{virus}_{protein}"
    aligns.setdefault(virus, {})
    aligns[virus].update({ref_seqid: aref})
    aligns[virus].update({query_seqid: aquery})

    # ensure lengths match up
    assert len(aquery) == len(aref) 
    # print("query and ref are the same length")
    
    # determine and attach coords
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
    loc = '.'.join('{}_{}'.format(r[0], r[1]) for r in ranges)
    coord_dict[virus][protein].update({"mapped_start": mapped_start, "mapped_stop": mapped_stop, "mapped_ranges": loc})
            
handle = open("/home/paula/surfaces_project/repos/surfaces/data/annotation_with_uniprot.csv", 'r')
outfile = open("/home/paula/surfaces_project/repos/surfaces/data/annotation_with_uniprot_step3_coords.csv", 'w')
reader = csv.DictReader(handle)
header = reader.fieldnames
header.extend(['mapped_start', 'mapped_stop', 'ranges'])

print(header)
writer = csv.DictWriter(outfile, fieldnames=header)
writer.writeheader()
for row in reader:
    virus = row['virus']
    protein = row['protein']
    # print(f"virus: protein") # debugging
    coord_data = coord_dict[virus][protein]
    row.update({"mapped_start": coord_data["mapped_start"], "mapped_stop": coord_data["mapped_stop"], "ranges": coord_data["mapped_ranges"]})

    writer.writerow(row)
handle.close()
outfile.close()


## SANITY CHECK

abbrev_dict = {
    "Chikungunya virus": "chikv",
    "dengue virus type 2": "dengue2",
    "Enterovirus A": "enteroA71/coxsacki",
    "Hepatitis C virus genotype 1": "HCV1a",
    "Human immunodeficiency virus 2": "HIV2",
    "Influenza A virus (A/Hong Kong/1073/99(H9N2))": "IAVH9N2",
    "Influenza B virus (B/Lee/1940)": "IBV",
    "Lyssavirus rabies":"lyssavirus",
    "Mamastrovirus 1": "Mastro1",
    "Measles morbillivirus": "measles",
    "Mumps orthorubulavirus": "mumps",
    "Borna disease virus 1": "orthobornavirus",
    "Enterovirus C":"polio",
    "Potato virus Y":"PotatoY",
    "Rhinovirus A": "RhinoA",
    "Rotavirus A": "RotaA",
    "Human respiratory syncytial virus A": "rsv",
    "Venezuelan equine encephalitis virus": "VEEV",
    "West Nile virus": "wnv",
    "Zika virus": "zika",
    "Human immunodeficiency virus 1": "HIV1A1",
    "Tick-borne encephalitis virus": "TBEV",
    "Foveavirus mali": "Foveavirus",
    "Tobacco mosaic virus": "Tobacco",
    "Potato virus X": "Potato",
    "Yellow fever virus": "YFV",
    "Hepatovirus A": "HepA"
}
# exceptions = ['Mamastrovirus 1'] # not sure what a good sub reference is for this

# identify the actual refseq coord for each virus/protein to check if ranges make sense
# need to account for polyprotein versus cds annotation
refseq_coords = {}
for accn in accns_set:
    record = SeqIO.read(f"/home/surfaces/1_rawdata/refseq/{accn}.gb", format="genbank")
    whole_nucseq = record.seq

    # get the virus name
    for feat in record.features:
        if feat.type == "source":
            virus = feat.qualifiers["organism"][0]
            abbrev_virus = abbrev_dict[virus] # get abbrev to align w file naming convention
            if abbrev_virus == "enteroA71/coxsacki": # this one appears to have used the same ref - identify more appopriate ref for enterov1 if possible
                abbrev_virus = "coxsacki"


    mat_peptides = [feat for feat in record.features if feat.type == 'mat_peptide']
    refseq_coords.setdefault(virus, {})
    if mat_peptides:
        # treat as mature peptides
        for mat_peptide in mat_peptides:
            q = mat_peptide.qualifiers

            parts = []
            for part in mat_peptide.location.parts:
                parts.append((part.start, part.end))
            loc = '.'.join('{}_{}'.format(p[0], p[1]) for p in parts)
            product = q.get('product', '')
            product = product[0]
            print(product)

            # get aligned segment seq for comp
            header = f"ref_{abbrev_virus}_{product}"
            segment_nucseq = str(mat_peptide.location.extract(record).seq)
            aquery, aref = mafft(segment_nucseq, whole_nucseq)


            aligns[abbrev_virus].update({header: aquery})

            # get exact coordinates for refseq proteins
            refseq_coords[virus].setdefault(product, {})
            refseq_coords[virus][product] = (loc, accn)
    
    else: # look at cds for proteins
        cds_annot = [feat for feat in record.features if feat.type == 'CDS']
        for cds in cds_annot:
            q = cds.qualifiers

            parts = []
            for part in cds.location.parts:
                parts.append((part.start, part.end))
            loc = '.'.join('{}_{}'.format(p[0], p[1]) for p in parts)

            product = q.get('product', '')
            product = product[0]
            print(product)

            header = f"ref_{abbrev_virus}_{product}"
            segment_nucseq = cds.location.extract(record).seq
            aquery, aref = mafft(segment_nucseq, whole_nucseq)
            aligns[abbrev_virus].update({header: aquery})

            # get exact coordinates for refseq proteins
            refseq_coords[virus].setdefault(product, {})
            refseq_coords[virus][product] = (loc, accn)

header = ["virus", "abbrev", "protein", "accession", "refseq_coords"]
with open("/home/paula/surfaces_project/repos/surfaces/data/refseq_step3_coords.csv", "w") as outfile_refseq:
    writer = csv.DictWriter(outfile_refseq, fieldnames=header)
    writer.writeheader()

    for virus, proteins in refseq_coords.items():
        for protein, data in proteins.items():
            writer.writerow({
                "virus": virus,
                "abbrev": abbrev_dict[virus],
                "protein": protein,
                "accession": data[1],
                "refseq_coords": data[0]
            })


# write out alignment - this should include whole genome, ref segments and aligned consesus fasta
for virus, seqs in aligns.items():
    seqout = open(f"/home/paula/surfaces_project/annotate/alignments/{virus}_step3_alignment.fasta", 'w')     
    for seqid, seq in seqs.items():
        seqout.write(f">{seqid}\n{seq}\n")
    seqout.close()