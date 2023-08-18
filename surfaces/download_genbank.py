import re
import os.path
from Bio import Entrez, SeqIO
from time import sleep
import csv

from mpi4py import MPI
the_world = MPI.COMM_WORLD
my_number = the_world.Get_rank()
total_number = the_world.Get_size()

Entrez.email = 'sbagher4@uwo.ca'

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
    return lst


def retrieve_record(accession_number):
    with Entrez.esearch(db="nucleotide", term=accession_number, retmax=1) as handle:
        response = Entrez.read(handle)
        gid = response['IdList'][0]

    with Entrez.efetch(db="nuccore", id=gid, rettype="gb", retmode="text") as handle:
        record = SeqIO.read(handle, "genbank")

    return record


def extract_coding_sequences(genbank_record, accession_number):
    coding_sequences = []
    for feature in genbank_record.features:
        if feature.type == "CDS":
            feature_dict = {}
            qualifiers = feature.qualifiers

            parts_list = []
            for part in feature.location.parts:
                parts_list.append((part.start, part.end))
                parts = ';'.join(f'{item[0].position}:{item[1].position}' for item in parts_list)

            feature_dict["parts"] = parts
            feature_dict["accession"] = accession_number
            feature_dict["protein_accession"] = qualifiers.get('protein_id', [''])[0]
            feature_dict["strand"] = feature.strand
            feature_dict["product"] = qualifiers.get('product', [''])[0]
            feature_dict["gene_sequence"] = feature.extract(genbank_record.seq)
            feature_dict["start"] = qualifiers['codon_start']

            coding_sequences.append(feature_dict)

    return coding_sequences


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


def retrieve_ref_cds_poly(accn_path, # File with all the ncbi accns
                          outdir):

    count = 0

    # list of accns
    accn_handle = open(accn_path, "r")
    accn_lst = accn_to_list(accn_handle)
    #sleep(1)

    for accn in accn_lst:
        count += 1
        if count % total_number != my_number:
            #sleep(1)
            continue

        dir_accn = os.path.join(outdir, accn)
        os.mkdir(dir_accn)

        record = retrieve_record(accn) # seqIO record
        #print(record)

        info = retrieve_poly_genes(record) # dict of genbank info
        #print(dict)

        cds_lst = extract_coding_sequences(accn, record) # #outfile name: '{}/{}_{}'.format(dir, accn, protein_accn)
        #print(cds_lst)

        # if its a polyprotien
        if info["poly"] == "TRUE":
            cds_d = cds_lst[0] # should only be one CDS listed for polyproteins
            # dict : {accn,cds_d["protein_accn"],cds_d["product"],cds_d["strand"],cds_d["parts_str"], cds_d["gene_sequence"]}

            # ref locations
            ref_start = cds_d["start"][0]

            # every mature peptide
            for p in info["protein_id"]:
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
                outfile.write(f'>"{accn},{prot_id},{info["product"]},{cds_d["strand"]},{location}"\n{all_cut}\n')
                print("poly " + accn + " " + prot_id)

	# not a polyprotein: d["poly"] = "FALSE"
        else:
            for cds_d in cds_lst:
		# write out all cds
                outfile = open('{}/{}_{}'.format(dir_accn, accn, cds_d["protein_accn"]), 'w')
                #>"NC_007545,YP_392488.1,nonstructural protein 2,1,42:981"
                outfile.write(f'>"{accn},{cds_d["protein_accn"]},{cds_d["product"]},{cds_d["strand"]},{cds_d["parts_str"]}"\n{cds_d["gene_sequence"]}\n')
                print("not poly " + accn + " " + cds_d["protein_accn"])


def parse_table(tbl):
    """
    Loop trough table to retrieve accession numbers as a list
    """
    handle = open(tbl)

    values = []
    for line in handle:
        temp = line.strip('\n')  # Remove new line
        values.append(temp)

    return(values)


def retrieve_CDS(record):
    """
    Analyze features in Genbank record to extract (1) the number of coding
    regions (CDS)
    """
    cds = [feat for feat in record.features if feat.type=='CDS']
    for cd in cds:
        q = cd.qualifiers
        parts = []
        for part in cd.location.parts:
            parts.append((part.start, part.end))
        locus = q.get('locus_tag', '')
        product = q.get('product', [''])
        aaseq = q.get('translation', [''])
        yield locus, product, cd.strand, parts, q['codon_start'], aaseq


def pase_filtered_neighbour_table(table):
    """
    #:param table: filtered_Neighbours.txt
    #:return: list of all the neighbour accessions
    """
    lst=[]
    with open (table, 'r') as f:
        title = next(f)
        header = next(f)
        reader=csv.reader(f,delimiter='\t')
        for line in reader:
            lst.append(line[1])

    return(lst)


def retrieve_annotations(record):
    an = record.annotations
    an_dict = {'Length': len(record.seq) , 'Topology':an['topology'],
                'Taxonomy':an['taxonomy'], 'Molecule': an['molecule_type'],
                'Proteins': len(record.features)
                }
    return an_dict


if __name__ == '__main__':
    accession_number = "KF853231"
    genbank_record = retrieve_record(accession_number)
    coding_sequences = extract_coding_sequences(genbank_record)
