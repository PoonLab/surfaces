description = """
Sample codon sites from an alignment at random with or without replacement.
Write the result to a new file.
"""

import argparse
import random
import sys
import os


def parse_fasta(handle):
    """ Parse headers and sequences from a FASTA file """
    headers = []
    sequences = []
    seq = ''

    for line in handle:
        if line.startswith(">"):
            if len(seq) > 0:
                # output preceding record
                headers.append(header)
                sequences.append(seq)
                seq = ''   # reset container
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
    columns = []
    for i in range(n_columns):
        l = step*i
        r = l + step
        columns.append([s[l:r] for s in seqs])
    return columns


def sample_codons(columns, k, replace=False):
    """
    Sample codon sites without replacement
    :param columns:  list of lists, each storing one codon site (column)
    :param k:  int, number of sites to sample
    :param replace:  bool, if True then sample with replacement
    :return:  list, strings produced by concatenating sampled columns
    """
    pop = range(len(columns))
    samp = random.choices(pop, k) if replace else random.sample(pop, k)
    res = [cod for cod in columns[samp[0]]]  # a list of strings

    for ci in samp[1:]:  # codon index
        for si in range(len(res)):  # sequence index
            res[si] += columns[ci][si]  # append codon to sequence
    return res
    

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
    
    for i in range(args.reps):
        res = sample_codons(columns, k=args.num, replace=args.replace)
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

