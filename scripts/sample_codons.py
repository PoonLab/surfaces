description = """
Sample codon sites from an alignment at random with or without replacement.
Write the result to a new file.
"""

from Bio import AlignIO, SeqIO
import argparse
import random
import sys
import os


def parse_fasta(handle):
    """ Parse headers and sequences from a FASTA file """
    header = None
    headers = []
    sequences = []
    seq = ''
    for line in handle:
        if line.startswith(">"):
            if len(seq) > 0:
                headers.append(header)
                sequences.append(seq)
                seq = ''   # reset containers
            header = line.strip('>\n')
        else:
            seq += line.strip()
    # handle last record
    headers.append(header)
    sequences.append(seq)

    return headers, sequences


def transpose_fasta(seqs, step=3):
    """ Convert list of sequences to list of columns """
    n_columns = len(seqs[0]) // step
    columns = [[] for i in range(n_columns)]
    for seq in seqs:
        for i in range(n_columns):
            columns[i].append(seq[(step*i):(step*(i+1))])
    return columns


def sample_codons(columns, k, reps=1, replace=False):
    """
    Sample codon sites without replacement
    :param columns:  list of lists, each storing one codon site (column)
    :return:  a new AlignIO.MultipleSeqAlignment object
    """
    for _ in range(reps):
        samp = random.choices(columns, k) if replace else random.sample(columns, k)
        res = samp[0]
        for block in samp[1:]:
            for i in range(len(res)):
                res[i] += block[i]  # append codon
        yield res
    

if __name__ == "__main__":
    # command line interface
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument("infile", type=argparse.FileType('r'), help="input, path to FASTA file")
    parser.add_argument("num", type=int, help="number of codon sites to sample")
    parser.add_argument("prefix", type=str, 
                        help="prefix (no extension) for path to write output.")

    parser.add_argument("--replace", '-r', action="store_true", 
                        help="option, set to sample with replacement")
    parser.add_argument("-f", "--format", type=str, default='fasta',
                        help="option, sequence file format (default: fasta)")
    parser.add_argument("--outfmt", type=str, default='fasta',
                        help="option, file format for output (default: fasta)")
    parser.add_argument("-n", "--reps", type=int, default=1, 
                        help="option, number of replicates (default: 1)")
    parser.add_argument("--overwrite", action="store_true",
                        help="set to overwrite existing output files.")
    parser.add_argument("--seed", type=int, default=None,
                        help="set random seed")

    args = parser.parse_args()
    assert args.reps > 0, "Error, -n must be greater than zero!"
    if args.seed:
        random.seed(args.seed)
    
    #msa = AlignIO.read(args.infile, args.format)
    #ln = msa.get_alignment_length() // 3
    headers, sequences = parse_fasta(args.infile)
    ln = len(sequences[0]) // 3

    if ln <= args.num:
        sys.stderr.write(f"Input alignment is already <= target length {args.num}\n")
        of = f"{args.prefix}_0.{args.outfmt}"
        with open(of, 'w') as handle:
            for h, s in zip(headers, sequences):
                handle.write(f">{h}\n{s}\n")
        sys.exit()

    columns = transpose_fasta(sequences, step=3)
    sampler = sample_codons(columns, k=args.num, reps=args.reps, replace=args.replace)
    
    for i, res in enumerate(sampler):
        of = f"{args.prefix}_{i}.{args.outfmt}"
        if os.path.exists(of) and not args.overwrite:
            sys.stderr.write(f"ERROR: file {of} already exists! "
                             "Set --overwrite to override.\n")
            sys.exit()
            
        with open(of, 'w') as handle:
            #SeqIO.write(res.sequences, handle, args.outfmt)
            #handle.write(format(res, args.outfmt))
            for h, s in zip(headers, res):
                handle.write(f">{h}\n{s}\n")

