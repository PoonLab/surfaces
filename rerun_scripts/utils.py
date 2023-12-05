import re
from Bio import Entrez, SeqIO

def parse_loc(item):
    """
    item : [949:1765](+)
    item = "join{[12312:12354](+), [12353:15131](+)}"
    item_normal = "[949:1765](+)"
    """
    matches = re.findall(r'\d+', item)
    if len(matches) == 2:
        start, stop = map(lambda x: int(x), matches)
        return start, stop


def get_accessions(handle):
    """
    input: handle to open file of accns separated by "\n"
    return: list of all the accessions
    ex. ['KF853231', 'DQ479956', 'JN704703']
    """
    accessions = []
    for line in handle:
        name = line.strip()
        accessions.append(name)
    return accessions


def retrieve_genbank(gid):
    """
    TODO
    """
    handle = Entrez.efetch(db='nuccore', id=gid, rettype='gb', retmode='text')
    gb = SeqIO.parse(handle, format='genbank')
    return next(gb)


def retrieve_record(accession_number):
    """
    TODO
    """
    with Entrez.esearch(db="nucleotide", term=f'{accession_number}[ACCN]', retmax=1) as handle:
        response = Entrez.read(handle)
        gid = response['IdList'][0]

    with Entrez.efetch(db="nuccore", id=gid, rettype="gb", retmode="text") as handle:
        record = SeqIO.read(handle, "genbank")

    return record


def retrieve_cds(accession_number, genbank_record):
    """
    input: accn = ncbi accn
    input: record = genbank file SeqIO.read
    input: dir = directory to write out cds
    return list of dictionary {accn},{protein_accn},{product},{strand},{parts_str}"\n{gene_sequence}
    output: [{'accn': 'NC_005300', 'protein_accn': 'NP_950235.1', 'product': 'glycoprotein precursor', 'strand': 1, 'parts_str': '92:5147', 'gene_sequence': Seq('ATGCATATATCATTAATGTATGCAATCCTTTGCCTACAGCTGTGTGGTCTGGGA...TAG')}]
    """
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

