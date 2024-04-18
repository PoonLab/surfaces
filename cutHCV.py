from Bio import SeqIO
import sys
import subprocess
import tempfile
import re
from io import StringIO


if len(sys.argv) != 2:
    print("Usage: python3 cutHCV.py <infile>")
    sys.exit()


codon_dict = {'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L',
                'TCT':'S', 'TCC':'S', 'TCA':'S', 'TCG':'S',
                'TAT':'Y', 'TAC':'Y', 'TAA':'*', 'TAG':'*',
                'TGT':'C', 'TGC':'C', 'TGA':'*', 'TGG':'W',
                'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L',
                'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P',
                'CAT':'H', 'CAC':'H', 'CAA':'Q', 'CAG':'Q',
                'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R',
                'ATT':'I', 'ATC':'I', 'ATA':'I', 'ATG':'M',
                'ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T',
                'AAT':'N', 'AAC':'N', 'AAA':'K', 'AAG':'K',
                'AGT':'S', 'AGC':'S', 'AGA':'R', 'AGG':'R',
                'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V',
                'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',
                'GAT':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E',
                'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G',
                '---':'-', 'XXX':'?'}

mixture_regex = re.compile('[WRKYSMBDHVN-]')

mixture_dict = {'W':'AT', 'R':'AG', 'K':'GT', 'Y':'CT', 'S':'CG',
                'M':'AC', 'V':'AGC', 'H':'ATC', 'D':'ATG',
                'B':'TGC', 'N':'ATGC', '-':'ATGC'}


def translate_nuc(seq, offset, resolve=False, return_list=False):
    """
    Translate nucleotide sequence into amino acid sequence.
            offset by X shifts sequence to the right by X bases
    Synonymous nucleotide mixtures are resolved to the corresponding residue.
    Nonsynonymous nucleotide mixtures are encoded with '?'
    """
    seq = '-' * offset + seq
    aa_list = []
    aa_seq = ''  # use to align against reference, for resolving indels

    # loop over codon sites in nucleotide sequence
    for codon_site in range(0, len(seq), 3):
        codon = seq[codon_site:codon_site + 3]
        if len(codon) < 3:
            break

        # note that we're willing to handle a single missing nucleotide as an ambiguity
        if codon.count('-') > 1 or '?' in codon:
            if codon == '---':  # don't bother to translate incomplete codons
                aa_seq += '-'
                aa_list.append(['-'])
            else:
                aa_seq += '?'
                aa_list.append(['?'])
            continue

        # look for nucleotide mixtures in codon, resolve to alternative codons if found
        num_mixtures = len(mixture_regex.findall(codon))
        if num_mixtures == 0:
            aa = codon_dict[codon]
            aa_seq += aa
            aa_list.append([aa])
        elif num_mixtures == 1:
            resolved_AAs = []
            for pos in range(3):
                if codon[pos] in mixture_dict.keys():
                    for r in mixture_dict[codon[pos]]:
                        rcodon = codon[0:pos] + r + codon[(pos + 1):]
                        if codon_dict[rcodon] not in resolved_AAs:
                            resolved_AAs.append(codon_dict[rcodon])

                        aa_list.append(resolved_AAs)

                        if len(resolved_AAs) > 1:
                            if resolve:
                                # for purposes of aligning AA sequences
                                # it is better to have one of the resolutions
                                # than a completely ambiguous '?'
                                aa_seq += resolved_AAs[0]
                            else:
                                aa_seq += '?'
                        else:
                            aa_seq += resolved_AAs[0]
                else:
                    aa_seq += '?'
                    aa_list.append(['?'])
    if return_list:
        return aa_list
    return aa_seq


def mafft(query, ref, binpath="mafft"):
    """
    Align translation of query sequence to a reference protein
    and use the resulting alignment to cut the corresponding region
    in the original nucleotide sequence of the query.

    :param query:  str, nucleotide sequence (polyprotein)
    :param ref:  str, reference sequence (amino acid)
    :param binpath:  str, path to MAFFT executable file
    :return:
    """
    tquery = translate_nuc(query, 0)  # aa translation

    handle = tempfile.NamedTemporaryFile(delete=False)
    s = f'>ref\n{ref}\n>query\n{tquery}\n'
    handle.write(s.encode('utf-8'))
    handle.close()

    # call MAFFT on temporary file
    stdout = subprocess.check_output([binpath, '--quiet', handle.name])
    stdout = stdout.decode('utf-8')
    result = list(SeqIO.parse(StringIO(stdout), "fasta"))
    aref = str(result[0].seq)
    aquery = str(result[1].seq)

    left = len(aref) - len(aref.lstrip('-'))
    right = len(aref) - len(aref.rstrip('-'))

    # cut corresponding interval from nucleotide sequence
    return query[(3*left):(3*right)]

handle = open(sys.argv[1])

# handle reference
ref = SeqIO.read(open("h77.gb"), "genbank")

proteins = {}
for feat in ref.features:
    if feat.type != 'mat_peptide':
        continue
    aaseq = str(feat.translate(ref).seq)
    product = feat.qualifiers['product'][0]
    proteins.update({product: aaseq})


records = SeqIO.parse(handle, 'fasta')
#aligned = dict([(p, {}) for p in proteins])
outfiles = {}
for p in proteins:
    pname = re.sub("[ /.,:]", "_", p)
    print(pname)
    outfile = open(f"{pname}.fasta", 'w')
    outfiles.update({p: outfile})

for record in records:
    query = str(record.seq)
    print(record.name)
    for protein, refseq in proteins.items():
        result = mafft(query=query, ref=refseq)
        if len(result) == 0:
            continue
        outfiles[protein].write(f">{record.description}\n{result}\n")

