from Bio import SeqIO
import argparse
import subprocess
import tempfile
from io import StringIO
import os
import sys


description = """
Use a reference Genbank record to classify CDS records by 
pairwise alignment.
"""


def parse_gb(gb_file):
    prots = {}
    record = SeqIO.read(gb_file, 'genbank')
    for feat in record.features:
        if feat.type == "CDS":
            prot = feat.qualifiers['translation'][0]
            product = feat.qualifiers['product'][0]
            prots.update({product: prot})
    return prots


def mafft(query, ref, binpath="mafft", match=1, mismatch=-1, gap=-3):
    """
    :param query:  str, sequence to align to reference
    :param ref:  str, reference sequence
    :param binpath:  str, path to MAFFT executable file
    :return:  int, alignment score
    """
    handle = tempfile.NamedTemporaryFile(delete=False)
    s = f'>ref\n{ref}\n>query\n{query}\n'
    handle.write(s.encode('utf-8'))
    handle.close()

    stdout = subprocess.check_output([binpath, '--quiet', handle.name])
    stdout = stdout.decode('utf-8')
    result = list(SeqIO.parse(StringIO(stdout), "fasta"))
    aref = str(result[0].seq)
    aquery = str(result[1].seq)
    
    ascore = 0
    for i, raa in enumerate(aref):
        if raa == '-':
            ascore += gap
            continue
        qaa = aquery[i]
        if qaa == '-':
            ascore += gap
        else:
            ascore += (mismatch if qaa != raa else match)            
    return ascore


if __name__ == "__main__":
    # command line interface
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("gb", type=argparse.FileType('r'), 
                        help="Genbank file for annotated reference genome")
    parser.add_argument("cds", type=argparse.FileType('r'), 
                        help="File containing CDS records")
    parser.add_argument('-f', '--format', type=str, default='fasta',
                        help="option, specify file format for CDS")
    parser.add_argument('--outdir', type=str, default='.',
                        help="option, dir for output files, e.g., data/")
    args = parser.parse_args()

    refs = parse_gb(args.gb)
    sys.stderr.write(f"Detected the following products: {','.join(list(refs.keys()))}")

    # create output files
    outfile = {}
    try:
        virus = os.path.basename(args.cds.name).split('_')[0]
    except:
        sys.stderr.write(f"Failed to parse virus name from {os.path.basename(args.cds.name)}")
        raise

    for key in refs:
        fn = f"{virus}_{key.split(' ')[0]}_step1.fasta" 
        outfile.update({key: open(os.path.join(args.outdir, fn), 'w')})

    for record in SeqIO.parse(args.cds, args.format):
        query = str(record.seq.translate())
        max_key = None
        max_score = -1e6
        for key, rseq in refs.items():
            ascore = mafft(query, rseq)
            if ascore > max_score:
                max_score = ascore
                max_key = key
        #print(record.description, max_key)
        SeqIO.write(record, outfile[max_key], 'fasta')

