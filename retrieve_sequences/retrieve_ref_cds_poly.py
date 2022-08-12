"""
handles polyproteins and non-polyproteins
takes in a text file with ncbi accns
output fasta files with all ref cds

FASTA HEADER: >"NC_010956,YP_001974451.1,core protein V,1,15696:16701"

FOR TESTING: 
accn_lst = ["NC_001526","NC_001451","NC_005300","NC_005222"]
dir = "/home/sareh/learning/test/cds/"
"""

import os
import re
import csv
from Bio import Entrez, SeqIO, SeqFeature
from csv import DictWriter
from time import sleep
from glob import glob
import argparse

Entrez.email = 'sbagher4@uwo.ca'

# mpi
from mpi4py import MPI
the_world = MPI.COMM_WORLD
my_number = the_world.Get_rank()
total_number = the_world.Get_size()

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--accn', type=str,
                       help='IN: file of accns (each line is an ncbi accn)')
    parser.add_argument('--dir', type=str,
                       help='OUT: directory to write out ref CDS')
    return parser.parse_args()


def parse_loc(item):
    """
    item : [949:1765](+)
    item = "join{[12312:12354](+), [12353:15131](+)}"
    item_normal = "[949:1765](+)"
    """
    loc = item.replace("(+)","").replace('[','').replace(']','')
    cord = loc.split(":")
    start = cord[0].strip(" ")
    stop = cord[1].strip(" ")
    return(start, stop)

def accn_to_list(handle):
    """
    input: handle to open file of accns separated by "\n"
    return: list of all the accessions
    ex.['KF853231', 'DQ479956', 'JN704703']
    """
    lst = []
    next(handle)
    for line in handle:
        name = line.split("\n")
        lst.append(name[0])
    return (lst)

def retrieve_gid(accn):
    """
    Find records according to genbank accession number
    """
    handle = Entrez.esearch(db='nucleotide', term=accn, retmax=1)
    response = Entrez.read(handle)
    if response['Count'] == '0':
        # retry query
        handle = Entrez.esearch(db='nucleotide', term=accn, retmax=1)
        response = Entrez.read(handle)
        if response['Count'] == '0':
            return None

    return response['IdList'][0]

def retrieve_genbank(gid):
    #handle = Entrez.efetch(db='nucleotide', rettype='gb', retmode='text', id=gid)
    handle  = Entrez.efetch(db="nuccore", id=gid, rettype="gb",retmode="text")
    gb = SeqIO.parse(handle, format='genbank')
    return next(gb)

def retrieve_record(accn):
    gid = retrieve_gid(accn)
    #sleep(1)
    if gid is None:
        print('Warning, failed to retrieve gid for {}'.format(accn))
        #continue
    record = retrieve_genbank(gid)
    #sleep(1) # don't spam, add more if it's a big file
    return(record)

def retrieve_cds(accn, record):
    """
    input: accn = ncbi accn
    input: record = genbank file SeqIO.read
    input: dir = directory to write out cds
    return list of dictionary {accn},{protein_accn},{product},{strand},{parts_str}"\n{gene_sequence}
    output: [{'accn': 'NC_005300', 'protein_accn': 'NP_950235.1', 'product': 'glycoprotein precursor', 'strand': 1, 'parts_str': '92:5147', 'gene_sequence': Seq('ATGCATATATCATTAATGTATGCAATCCTTTGCCTACAGCTGTGTGGTCTGGGA...TAG')}]
    """
    # retrieve each CDS
    cds = [feat for feat in record.features if feat.type == "CDS"]  # [SeqFeature(FeatureLocation(ExactPosition(154), ExactPosition(544), strand=1), type='CDS'), SeqFeature(FeatureLocation(ExactPosition(570), ExactPosition(897), strand=-1), type='CDS'), SeqFeature(FeatureLocation(ExactPosition(1115), ExactPosition(1832), strand=1), type='CDS')

    # dict of cds
    list_cds_dict = []

    for cd in cds:
        d = {}
        q = cd.qualifiers  # OrderedDict([('gene', ['ORF 0']), ('codon_start', ['1']), ('product', ['ORF 0']), ('protein_id', ['AHJ09141.1']), ('translation', ['MATVHYSRRPGTPPVTLTSSPGMDDVATPIPYLPTYAEAVADAPPPYRSRESLVFSPPLFPHVENGTTQQSYDCLDCAYDGIHRLQLAFLRIRKCCVPAFLILFGILTLTAVVVAIVAVFPEEPPNSTT'])])
        proteinID = q.get('protein_id', [''])  # ['AHJ09212.1']
        protein_accn = proteinID[0]  # AHJ09141.1
        product = q.get('product', [''])[0]  # ORF 0
        locus = q.get('locus_tag', '')  # empty
        strand = cd.strand  # 1 or -1
        start = q['codon_start']  # ['1']
        seq = record.seq  # full nucleotide seqeunce
        name = record.annotations['organism']
        gene_sequence = cd.extract(record.seq)

        parts = []  # [(ExactPosition(154), ExactPosition(544))]
        for part in cd.location.parts:  # [154:544](+)
            parts.append((part.start, part.end))
            parts_str = ';'.join('{}:{}'.format(p[0].position, p[1].position) for p in parts)

        d["accn"] = accn
        d["protein_accn"] = protein_accn
        d["product"] = product
        d["strand"] = strand
        d["parts_str"] = parts_str
        d["gene_sequence"] = gene_sequence
        d["start"] = start
        #outfile = open('{}/{}_{}'.format(dir, accn, protein_accn), 'w')
        #outfile.write(f'>"{accn},{protein_accn},{product},{strand},{parts_str}"\n{gene_sequence}\n')

        list_cds_dict.append(d)

        ### NEED to loop through and return all ###
        sleep(1)

    return(list_cds_dict)

def retrieve_poly_genes(record):
    """
    polyproteins only have one CDS listed
    so need to get the mature peptides (each a gene)
    retrun a dictionary [id:ncbi_accn, poly:boolean, seg:boolean, protien_id:[], CDS:[]]
    """
    d = {} # a new dicts for every accn
    d["id"] = record.id
    d["poly"] = "FALSE"
    d["seg"] = "FALSE"
    d["protein_id"] = []

    for feature in record.features:
        # check if segmented
        if feature.type == "source" and 'segment' in feature.qualifiers:
            d["seg"] = "TRUE"

        # check if its a polyprotein
        if feature.type == 'mat_peptide':
            d["poly"] = "TRUE"

            try:
                loc = feature.location

                #protein ID
                prot_accn = feature.qualifiers['protein_id'][0]
                d["protein_id"].append((prot_accn,str(loc)))

                #protein name
                product = feature.qualifiers['product'][0]
                d["product"] = product

                """
                #csv file only for polyprotein
                line = "{},{},{}\n".format(record.id, prot_accn, loc)
                loc_handle = open(args.loc, 'a')
                loc_handle.write(line)
                """

            except:
                # skip viruses with missing record ID
                print("issue: " + record.id)
    return(d)


def main():

    args = parse_args()
    count = 0

    accn_path = args.accn # File with all the ncbi accns

    # list of accns
    accn_handle = open(accn_path, "r")
    accn_lst = accn_to_list(accn_handle)
    #sleep(1)

    for accn in accn_lst:
        count += 1
        if count % total_number != my_number:
            #sleep(1)
            continue

        dir_accn = os.path.join(args.dir,accn)
        os.mkdir(dir_accn)

        record = retrieve_record(accn) # seqIO record
        #print(record)

        dict = retrieve_poly_genes(record) # dict of genbank info
        #print(dict)

        cds_lst = retrieve_cds(accn, record) # #outfile name: '{}/{}_{}'.format(dir, accn, protein_accn)
        #print(cds_lst)

        # if its a polyprotien
        if dict["poly"] == "TRUE":
            cds_d = cds_lst[0] # should only be one CDS listed for polyproteins
            # dict : {accn,cds_d["protein_accn"],cds_d["product"],cds_d["strand"],cds_d["parts_str"], cds_d["gene_sequence"]}

            # ref locations
            ref_start = cds_d["start"][0]

            # every mature peptide
            for p in dict["protein_id"]:
                prot_id = p[0] #NP_740469.1
                accn_full = accn
                accn = re.sub("\.[1-9]","",accn_full) #NC_002058: remove version
                location = p[1] #[949:1765](+)

                location = location.replace("join","").replace("{","").replace("}","")
                loc_lst = location.split(",")

                #ISSUE: join{[12312:12354](+), [12353:15131](+)}
                all_cut = ""
                for loc in loc_lst:
                    start,stop = parse_loc(loc) #poly start/stop

                    # cutting it out
                    cut_start = int(start) - int(ref_start)
                    cut_stop = int(stop) - int(ref_start) + 1
                    cut = record.seq[cut_start:cut_stop]
                    all_cut += str(cut) #add the string

		# fasta output
                outfile = open('{}/{}_{}'.format(dir_accn, accn, prot_id), 'w')
                #>"NC_007545,YP_392488.1,nonstructural protein 2,1,42:981"
                outfile.write(f'>"{accn},{prot_id},{dict["product"]},{cds_d["strand"]},{location}"\n{all_cut}\n')
                print("poly " + accn + " " + prot_id)

	# not a polyprotein: d["poly"] = "FALSE"
        else:
            for cds_d in cds_lst:
		# write out all cds
                outfile = open('{}/{}_{}'.format(dir_accn, accn, cds_d["protein_accn"]), 'w')
                #>"NC_007545,YP_392488.1,nonstructural protein 2,1,42:981"
                outfile.write(f'>"{accn},{cds_d["protein_accn"]},{cds_d["product"]},{cds_d["strand"]},{cds_d["parts_str"]}"\n{cds_d["gene_sequence"]}\n')
                print("not poly " + accn + " " + cds_d["protein_accn"])

if __name__ == "__main__":
    main()
