"""
For a virus genome that encodes a single polyprotein, we want to partition 
this sequence up into individual gene products (mature peptides).
This script will probably not work for viruses with gene products that 
result from spliced transcripts (joined from different intervals of the 
genome sequence).  For example, HIV-1 rev.
"""

from Bio import SeqIO
import sys
import subprocess
import tempfile
import re
import argparse
from io import StringIO
import pprint

# codon to amino acid conversion
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

def get_args(parser):
    # Required arguments
    parser.add_argument(
        'ref', type=argparse.FileType('r'),
        help='Genbank file with mature peptide'
    )
    parser.add_argument(
        'polyprots',
        help='fasta file with polyprotein encoding CDSs'
    )
    # Optional
    parser.add_argument(
        '--label', default=None, 
        help='<Optional> Path to output files with CDS and \
                amino acid sequences of mature peptides'
    )
    parser.add_argument(
        '--rf', default=1, type=int ,
        help='<Optional> reading frames to evaluate. Can be 1, 2 or 3'
    )
    parser.add_argument(
        '--gt', default=0.5, type=float,
        help='<Optional> gap threshold: maximum percentage of internal gaps\
              allowed in the refseq when aligning amino acids'
    )
    return parser.parse_args()


def translate_nuc(seq, offset, resolve=False, return_list=False):
    seq = '-' * offset + seq.upper()
    aa_seq = ''
    for codon_site in range(0, len(seq), 3):
        codon = seq[codon_site:(codon_site+3)]
        if len(codon) < 3:
            break
        check = sum([codon.count(nt) for nt in 'ACGT'])
        aa_seq += codon_dict[codon] if check == 3 else 'X'
    return aa_seq


def mafft(query, ref, binpath='mafft', gap_open=10.0):
    """ Wrapper for MAFFT """
    handle = tempfile.NamedTemporaryFile(delete=False)
    s = f'>ref\n{ref}\n>query\n{query}\n'
    handle.write(s.encode('utf-8'))
    handle.close()

    # call MAFFT on temporary file
    stdout = subprocess.check_output(
        [binpath, '--op', str(gap_open), '--quiet', handle.name])
    stdout = stdout.decode('utf-8')
    result = list(SeqIO.parse(StringIO(stdout), "fasta"))
    aref = str(result[0].seq)
    aquery = str(result[1].seq)
    return aquery, aref


def align_score(aquery, aref, match=1, mismatch=-1, gap=-2):
    """ A simple pairwise alignment score """
    if len(aquery) != len(aref):
        return None
    score = 0
    for i, raa in enumerate(aref):
        qaa = aquery[i]
        if raa == '-' or qaa == '-':
            score += gap
        else:
            score += match if raa == qaa else mismatch
    score /= len(aref)
    return score


def apply_align(nucseq, aquery, left, right):
    """ Apply pairwise amino acid alignment to nucleotides """
    new_seq = ''
    qpos = 3*left
    for i in range(left, right):
        if aquery[i] == '-':
            continue  # skip insertion in reference
        new_seq += nucseq[qpos:(qpos+3)]
        qpos += 3
    return new_seq


def extract(query, ref, binpath="mafft", rf=1, gap_threshold = 0.5, gap_open = 10.0):
    """
    Align translation of query sequence to a reference protein
    and use the resulting alignment to cut the corresponding region
    in the original nucleotide sequence of the query.

    :param query: str, nucleotide sequence (polyprotein)
    :param ref: str, reference sequence (amino acid)
    :param binpath: str, path to MAFFT executable file
    :param rf: int, number of reading frames for protein translation (i.e. 1, 2 or 3)
    :param gap_percent: int, number of reading frames for protein translation
    :return:
    """
    # translate input nucleotide sequence into amino acids
    aln_scores = {}
    align_cache = {}
    for rf in range(0, rf):
        tquery = translate_nuc(query, rf)
        aquery, aref = mafft(tquery, ref)

        # Extract coordinates in reference
        left = len(aref) - len(aref.lstrip('-'))
        right = len(aref.rstrip('-'))
        refseq = aref[left:right]
        qseq = aquery[left:right]

        # score = No of matches - number of gaps within refseq
        int_gaps = refseq.count('-')  # internal gaps
        ascore = align_score(qseq, refseq)
        
        # Store nucleotide sequences bellow gap threshold
        if (int_gaps/len(refseq)) > gap_threshold:
            sys.stderr.write(f"Bad alignment:\n{refseq}\n{aquery[left:right]}\n")
            sys.exit()
        aln_scores[rf] = ascore
        aln_nuc = apply_align(query, aquery, left, right)
        align_cache[rf] = {'nuc': aln_nuc, 'aquery': qseq, 'aref': refseq}

    if aln_scores:
        # Find the best alignment 
        best_score = max(aln_scores.values())
        best_rf = max(aln_scores, key=aln_scores.get)
        max_nt_seq = align_cache[best_rf]['nuc']
        if best_score < 0:
            sys.stderr.write(f"Potential bad alignment with score {best_score}:\n" \
                f"{align_cache[best_rf]['aref']}\n{align_cache[best_rf]['aquery']}\n")
    else:
        max_nt_seq = ''

    return max_nt_seq


def get_reference_proteins(ref_genome):
    """
    from a Genbank file a virus with map_peptide annotation, 
    extract the proteins
    :param ref_genome: SeqIO object
    :return proteins: dict, keyed by protein time, 
                            values are amino acid sequences
    """
    proteins = {}
    for feat in ref_genome.features:
        if feat.type != 'mat_peptide':
            continue
        aaseq = str(feat.translate(ref).seq)
        product = feat.qualifiers['product'][0]
        proteins.update({product: aaseq})
    return (proteins)


if __name__=="__main__":
    parser = argparse.ArgumentParser(
        description = "Extract mature peptides by mapping polyprotein to reference genome"
    )
    args = get_args(parser)

    polyprots = args.polyprots
    ref = SeqIO.read(args.ref, "genbank")
    label = args.label if args.label else ''

    # extract mature peptide sequences from reference
    proteins = get_reference_proteins(ref)
    proteins = {"1D": proteins["1D"]}  # testing

    # prepare separate output FASTA files for different genes/proteins
    outfiles = {}
    for p in proteins:
        pname = re.sub("[/.,:]", "_", p)
        pname = pname.split()[0]  #re.sub("_mature_peptide", "", pname)
        print(pname)  # debugging
        outfile = open(f"{label}_{pname}_step1.fasta", 'w')
        outfiles.update({p: outfile})
    
    # load the input FASTA file
    records = SeqIO.parse(polyprots, 'fasta')
    not_found = {}
    for record in records:
        query = str(record.seq)  # CDS file contains entire polyprotein sequences (in-frame!)
        name_parts = record.name.split("-")
        prot_name = name_parts[0]

        # Search each protein in reference genome
        for protein, refseq in proteins.items():
            name_parts[2] = re.sub("[ /.,:]", "_", protein)
            # Get the part of nucleotide sequence
            result = extract(query=query, ref=refseq, rf=args.rf, gap_threshold=args.gt)
            if len(result) == 0:
                # print(f"\n>>>>> no result for protein {protein} in {prot_name}<<<<<\n")
                if prot_name not in not_found.keys():
                    not_found[prot_name] = []
                not_found[prot_name].append(protein)
                continue
            header = "-".join(name_parts)
            outfiles[protein].write(f">{record.name}\n{result}\n")
    
    print("Not found proteins:")
    pprint.pprint(not_found)
