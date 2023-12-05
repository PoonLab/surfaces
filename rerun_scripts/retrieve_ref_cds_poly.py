import os
import re
from Bio import Entrez
import argparse

from utils import *

Entrez.email = 'sbagher4@uwo.ca'

# mpi
from mpi4py import MPI
the_world = MPI.COMM_WORLD
my_number = the_world.Get_rank()
total_number = the_world.Get_size()

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--infile', type=str,
                        help='IN: file of accns (each line is an ncbi accn)')
    parser.add_argument('--outdir', type=str,
                        help='OUT: directory to write out ref CDS')
    return parser.parse_args()


def retrieve_poly_genes(record):
    """
    polyproteins only have one CDS listed
    so need to get the mature peptides (each a gene)
    retrun a dictionary [id:ncbi_accn, poly:boolean, seg:boolean, protien_id:[], CDS:[]]
    """
    poly_genes = {} # a new dict for every accn
    poly_genes["id"] = record.id
    poly_genes["poly"] = False
    poly_genes["seg"] = False
    poly_genes["protein_id"] = []
    poly_genes["product"] = []

    for feature in record.features:
        poly_genes["seg"] = feature.type == "source" and 'segment' in feature.qualifiers
        poly_genes["poly"] = feature.type == 'mat_peptide'

        if poly_genes["poly"]:
            try:
                loc = feature.location #loc 
                prot_accn = feature.qualifiers['protein_id'][0] # protein ID
                product = feature.qualifiers['product'][0] # protein name
                poly_genes["protein_id"].append((prot_accn,str(loc),product))
            except:
                print("issue: " + record.id)

    return poly_genes


def retrieve_ref_cds(accn_path, accn_outdir):
    mpi_count = 0 # for mpi run
    ref_accn_count = 0
    poly_genes = []
    
    accn_handle = open(accn_path, "r")
    accessions = get_accessions(accn_handle)

    for accession in accessions:
        ref_accn_count += 1
        mpi_count += 1
        if mpi_count % total_number != my_number:
            continue

        outdir = os.path.join(accn_outdir, accession)

        if os.path.exists(outdir):
            continue

        os.mkdir(outdir)
        record = retrieve_record(accession)
        genbank_info = retrieve_poly_genes(record)
        cds_lst = retrieve_cds(accession, record)

        if genbank_info["poly"]:
            cds_d = cds_lst[0] # should only be one CDS listed for polyproteins
            ref_start = cds_d["start"][0]

            # every mature peptide
            for protein in genbank_info["protein_id"]:
                prot_id, location, prot_name = protein
                poly_genes.append(prot_id)
                accession = re.sub("\.[1-9]", "", accession)

                location = location.replace("join","").replace("{","").replace("}","")
                loc_lst = location.split(",")

                all_cut = ""
                for loc in loc_lst:
                    start, stop = parse_loc(loc)
                    cut_start = int(start) - int(ref_start)
                    cut_stop = int(stop) - int(ref_start) + 1
                    cut = record.seq[cut_start:cut_stop]
                    all_cut += str(cut)

                outfile = open(f'{outdir}/{accession}_{prot_id}', 'w')
                outfile.write(f'>"{accession},{prot_id},{prot_name},{cds_d["strand"]},{location}"\n{all_cut}\n')

        else:
            for cds_d in cds_lst:
                outfile = open('{}/{}_{}'.format(outdir, accession, cds_d["protein_accession"]), 'w')
                outfile.write(f'>"{accession},{cds_d["protein_accession"]},{cds_d["product"]},{cds_d["strand"]},{cds_d["parts"]}"\n{cds_d["gene_sequence"]}\n')

    print(poly_genes)


if __name__ == "__main__":
    args = parse_args()
    retrieve_ref_cds(args.infile, args.outdir)
