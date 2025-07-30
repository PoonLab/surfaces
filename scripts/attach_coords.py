from Bio import SeqIO
import csv
import os
from glob import glob
from Bioplus.consensus import conseq
from Bioplus.pair_align import mafft
import argparse
import sys

description = """
Use reference Genbank record to map consensus sequence by
pairwise alignment
"""

def parse_annot(reader):
    # parse annotation CSV
    accns = {}
    accns_set = set()
    for row in reader:
        virus = row['virus']
        if virus not in accns:
            accns.update({virus: {}})
        accn = row['refseq'].split('.')[0]
        accns[virus].update({row['protein']: accn})
        accns_set.add(accn)
    return accns, accns_set

def retrieve_refseq(gb_file):
    # retrieve sequences from Genbank files
    record = SeqIO.read(gb_file, 'genbank')
    accn = os.path.basename(gb_file).split('.')[0]
        
    for feat in record.features:
        if feat.type == "source":
            seg_no = feat.qualifiers.get("segment", [None])[0]
            if seg_no is None:
                break
    return accn, record.seq, seg_no

def map_consensus(accn, virus, protein):
    # retrieve reference
    refseq = refseqs.get(accn, None)
    if refseq is None:
        sys.stderr.write(f"WARNING: failed to retrieve reference for {accn}\n")
        return
    
    # generate consensus sequence
    aln = SeqIO.parse(fn, 'fasta')
    cseq = conseq(aln, thresh=0.5, is_nuc=True) # Bioplus
        
    # align to corresponding reference
    aquery, aref = mafft(cseq, refseq)

    if virus in segmented:
        ref_seqid = f"ref_{virus}_full"
        query_seqid = f"{virus}_{protein}"
        aligns.setdefault(virus, {}).setdefault(accn, {})
        aligns[virus][accn].update({ref_seqid: aref})
        aligns[virus][accn].update({query_seqid: aquery})

    else:
        ref_seqid = f"ref_{virus}_full"
        query_seqid = f"{virus}_{protein}"
        aligns.setdefault(virus, {})
        aligns[virus].update({ref_seqid: aref})
        aligns[virus].update({query_seqid: aquery})

    # ensure lengths match up
    assert len(aquery) == len(aref) 
    
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
    return mapped_start, mapped_stop, loc

def map_reference(accn, refseq_coords, abbrev_dict, gb_path, redo, verbose=False):
    record = SeqIO.read(f"{gb_path}{accn}.gb", format="genbank")
    whole_nucseq = record.seq

    # get the virus name
    for feat in record.features:
        if feat.type == "source":
            virus = feat.qualifiers["organism"][0]
            abbrev_virus = abbrev_dict[virus] # get abbrev to align w file naming convention
            # this one appears to have used the same ref - find more appopriate ref for enterov1
            if abbrev_virus == "enteroA71/coxsacki": 
                abbrev_virus = "coxsacki"
            if redo:
                if abbrev_virus not in redo: #DEBUGGING
                    continue
    
    if verbose:
        print(f"mapping {abbrev_virus} to reference")

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

            if type(product) == list:
                if len(product) == 1:
                    product = product[0]
                else:
                    sys.stderr.write(
                        f"failed to map consensus for {virus}-{protein} to reference genome")
                    sys.exit()

            # get aligned segment seq for comp
            header = f"ref_{abbrev_virus}_{product}"
            segment_nucseq = str(mat_peptide.location.extract(record).seq)
            aquery, aref = mafft(segment_nucseq, whole_nucseq)

            if virus in segmented:
                aligns[abbrev_virus][accn].update({header: aquery})
            else:
                aligns[abbrev_virus].update({header: aquery})

            
            # get exact coordinates for refseq proteins
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
            if product == '':
                product = q.get('note', '')
            product = product[0]
            product = product.replace(' ', '_')

            header = f"ref_{abbrev_virus}_{product}"
            segment_nucseq = cds.location.extract(record).seq
            aquery, aref = mafft(segment_nucseq, whole_nucseq)

            if abbrev_virus in segmented:
                aligns[abbrev_virus][accn].update({header: aquery})

            else:
                aligns[abbrev_virus].update({header: aquery})

            # get exact coordinates for refseq proteins
            refseq_coords[virus][product] = (loc, accn)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("-gb", "--gbpath", type=str,
                        default=r"/home/surfaces/1_rawdata/refseq/",
                        help="dir for reference genbank files")
    parser.add_argument("-a", "--annot", type=str, 
                        default="../data/annotation_with_uniprot.csv",
                        help="annotation csv")
    parser.add_argument("-f", "--fasta", type=str,
                        default="/home/surfaces/3_codaln",
                        help="path to fasta files for consensus seq")
    parser.add_argument("-o", "--outpath", type=str,
                        default="/home/paula/surfaces_project/annotate/test",
                        help="outfile path")
    parser.add_argument("--redo", action="store_true", 
                        default=False,help="pass list of viruses to redo")
    parser.add_argument("--verbose", action="store_true",
                        default=False, help="increase output verbosity")
    
    args = parser.parse_args()
    gb_path = args.gbpath
    fasta_path = args.fasta
    outpath = args.outpath
    verbose = args.verbose

    # parse annotation csv
    handle = open(args.annot, 'r')
    reader = csv.DictReader(handle)
    accns, accns_set = parse_annot(reader)
    handle.close()

    # retrieve reference sequence and segment info (if applicable) from Genbank files
    segments = {}
    refseqs = {}
    
    for fn in glob(f"{gb_path}*.gb"):
        accn, refseq, segment  = retrieve_refseq(fn)
        segments.update({accn: segment})
        refseqs.update({accn: refseq})

    # process fasta files
    types = ("*step3.fasta", "*step3.cds.fa", "measles*step2.fasta") # measles only has step 2 fasta
    files_3 = []
    for t in types:
        files_3.extend(glob(f"{fasta_path}/{t}"))

    redo = args.redo
    segmented = []
    if redo:
        for virus in segmented: # only run ones if in redo
            if virus in redo:
                segmented.append(virus)
    else:
        segmented = ["RotaA", "IBV", "IAVH9N2"]

    coord_dict = {}
    aligns = {}
    for fn in files_3:
        if args.verbose:
            print(f"processing fasta: {os.path.basename(fn)}")   

        tokens = os.path.basename(fn).split('_')
        virus = tokens[0]
        if redo:
            if virus not in redo:
                continue
            else:
                if args.verbose:
                    print(f"rerunning: {virus}")

        if virus not in accns:
            sys.stderr.write(f"ERROR: Could not find virus {virus} in annotations\n")
            sys.exit()
        
        protein = ' '.join(tokens[1:-1])
        # print(f"{virus}-{protein}")
        if protein not in accns[virus]:
            if ".fasta" in protein or ".fa" in protein: # handling misnamed files 
                protein = protein.split('.')[0]
            else:
                sys.stderr.write(
                    f"ERROR: Could not find {virus} protein '{protein}' in annotations\n")
                sys.exit()

        
        accn = accns[virus][protein]
        coord_dict.setdefault(virus, {})
        coord_dict[virus].setdefault(protein, {})

        # map consensus to reference genome
        start, stop, loc = map_consensus(accn, virus, protein)
        coord_dict[virus][protein].update({"mapped_start": start,"mapped_stop": stop,
                                           "mapped_ranges": loc})

        if not coord_dict:
            sys.stderr.write(
                    f"ERROR: failed to map consensus for {virus}-{protein} to reference genome\n")
            sys.exit()

    handle = open(args.annot, 'r')
    outfile = open(f"{outpath}/annotation_with_uniprot_step3_coords.csv", 'w')
    reader = csv.DictReader(handle)
    header = reader.fieldnames
    header.extend(['mapped_start', 'mapped_stop', 'ranges'])

    writer = csv.DictWriter(outfile, fieldnames=header)
    writer.writeheader()
    for row in reader:
        virus = row['virus']
        if redo:
            if virus not in redo: # DEBUGGING 
                continue
        protein = row['protein']
        coord_data = coord_dict[virus][protein]
        row.update({"mapped_start": coord_data["mapped_start"],
                    "mapped_stop": coord_data["mapped_stop"], "ranges": coord_data["mapped_ranges"]})

        writer.writerow(row)
    handle.close()
    outfile.close()


    ## SANITY CHECK - dict w abbrev
    abbrev_dict = {
        "Chikungunya virus": "chikv",
        "dengue virus type 2": "dengue2",
        "Enterovirus A": "enteroA71/coxsacki",
        "Hepatitis C virus genotype 1": "HCV1a",
        "Human immunodeficiency virus 2": "HIV2",
        "Influenza A virus (A/Hong Kong/1073/99(H9N2))": "IAVH9N2",
        "Influenza B virus (B/Lee/1940)": "IBV",
        "Lyssavirus rabies":"lyssavirus",
        # "Mamastrovirus 1": "Mastro1", # old ref
        "Human astrovirus 1": "Mastro1",
        "Measles morbillivirus": "measles",
        "Mumps orthorubulavirus": "mumps",
        "Borna disease virus 1": "orthobornavirus",
        "Enterovirus C":"polio",
        "Potato virus Y":"PotatoY",
        # "Rhinovirus A": "RhinoA", # old ref
        "rhinovirus A1": "RhinoA",
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

    refseq_coords = {}
    for accn in accns_set:
        map_reference(accn, refseq_coords, abbrev_dict, gb_path=gb_path, redo=redo, verbose=verbose)

    if not refseq_coords:
        sys.stderr.write(f"ERROR: reference mapping failed\n")
        sys.exit()

    header = ["virus", "abbrev", "protein", "accession", "refseq_coords"]
    with open(f"{outpath}/refseq_step3_coords.csv", 'w') as outfile_refseq:
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

    # write out alignment - include whole genome, ref segments and aligned consesus fasta
    for virus, seqs in aligns.items():
        if virus in segmented:
            for accn, seq in seqs.items():
                seg_no = segments[accn]
                seqout = open(f"{outpath}/{virus}_segment{seg_no}_step3_alignment.fasta", 'w')
                for seqid, nuc in seq.items():
                    seqout.write(f">{seqid}\n{nuc}\n")
                seqout.close() 
        else:
            seqout = open(f"{outpath}/{virus}_step3_alignment.fasta", 'w')     
            for seqid, seq in seqs.items():
                seqout.write(f">{seqid}\n{seq}\n")
            seqout.close()