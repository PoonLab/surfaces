description = """
Sample codon sites from an alignment at random with or without replacement.
Write the result to a new file.
"""

from Bio import AlignIO, SeqIO
import argparse
import random
import sys
import os
    

def sample_codons(msa, k, reps=1, replace=False):
    """
    Sample codon sites without replacement
    :param msa:  object of class AlignIO.MultipleSeqAlignment
    :return:  a new AlignIO.MultipleSeqAlignment object
    """
    ln = msa.get_alignment_length()
    if ln % 3 != 0:
        sys.stderr.write(f"Error in sample_codons(): alignment length ({ln}) is not "
                         "divisible by 3!\n")
        sys.exit()
    
    # partition alignment into codon sites
    columns = [msa.alignment[:, i:(i+3)] for i in range(0, ln, 3)]
    
    for _ in range(reps):
        samp = random.choices(columns, k) if replace else random.sample(columns, k)
        res = None
        for block in samp:
            if res is None:
                res = block
                continue
            res += block
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
    
    msa = AlignIO.read(args.infile, args.format)
    ln = msa.get_alignment_length() // 3
    if ln <= args.num:
        sys.stderr.write(f"Input alignment is already <= target length {args.num}\n")
        of = f"{args.prefix}_0.{args.outfmt}"
        with open(of, 'w') as handle:
            handle.write(format(msa, args.outfmt))
        sys.exit()
    
    sampler = sample_codons(msa, k=args.num, reps=args.reps, replace=args.replace)
    
    for i, res in enumerate(sampler):
        of = f"{args.prefix}_{i}.{args.outfmt}"
        if os.path.exists(of) and not args.overwrite:
            sys.stderr.write(f"ERROR: file {of} already exists! "
                             "Set --overwrite to override.\n")
            sys.exit()
            
        with open(of, 'w') as handle:
            #SeqIO.write(res.sequences, handle, args.outfmt)
            handle.write(format(res, args.outfmt))
