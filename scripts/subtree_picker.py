from Bio import SeqIO, Phylo
import argparse
import sys
from bisect import bisect

description = """
Down-sample sequences in a FASTA file based on a user-supplied tree.
Check that all labels match between these input files.
If --target and --branch are not provided, then output an alignment 
containing only sequences that match labels in the tree.
Otherwise, select the smallest subtree containing the --branch
label and related tips until the total length first exceeds --target.
"""

parser = argparse.ArgumentParser(description=description)
parser.add_argument("fasta", type=argparse.FileType('r'), 
                    help="input, path to FASTA file")
parser.add_argument("tree", type=argparse.FileType('r'),
                    help="input, path to file containing Newick tree")
parser.add_argument("--target", type=float, 
                    help="option, target length of subtree.")
parser.add_argument("--branch", type=str, 
                    help="option, label to situate subtree from")
parser.add_argument("--outfile", type=argparse.FileType('w'), default=sys.stdout,
                    help="option, path to write down-sampled FASTA (default: stdout)")
parser.add_argument("--midpoint", action="store_true", 
                    help="option, set to re-root the tree at its midpoint")
parser.add_argument("--csvfile", type=argparse.FileType('w'),
                    help="option, path to write CSV file of identical sequence labels")
args = parser.parse_args()

# import data from FASTA file
records = SeqIO.parse(args.fasta, "fasta")
fasta = {}
for record in records:
    fasta.update({record.description: record.seq})

# import tree and check labels against FASTA
tree = Phylo.read(args.tree, 'newick')
if args.midpoint:
    tree.root_at_midpoint()

branch = None
for tip in tree.get_terminals():
    if tip.name not in fasta:
        sys.stderr(f"ERROR: Could not locate {tip.name} in FASTA!\n")
        sys.exit(1)
    if args.branch and args.branch in tip.name:
        if branch is not None:
            sys.stderr.write("ERROR: found multiple hits in FASTA for --branch\n")
            sys.exit()
        else:
            branch = tip  # store Clade object

if args.branch and branch is None:
    sys.stderr.write(f"ERROR: --branch {args.branch} is not in tree!\n")
    sys.exit()


# by default, export only sequences that appear in the tree
if args.target is None or args.branch is None:
    sys.stderr.write("No target set by user, writing all tips in tree to FASTA by default...\n")
    for tip in tree.get_terminals():
        args.outfile.write(f">{tip.name}\n{fasta[tip.name]}\n")
    args.outfile.close()
    sys.exit()


assert args.target > 0, "--target must be greater than zero!"
assert args.target < tree.total_branch_length(), "--target exceeds length of entire tree!"

# otherwise traverse from user-specified tip toward root
nodes = tree.get_path(branch)
nodes.reverse()  # start with tip

# get lengths of all subtrees rooted at each node on path
subtrees = [tree.from_clade(node) for node in nodes]
tlens = [st.total_branch_length() for st in subtrees]

# find subtree with length closest to target
delta = [abs(tl-args.target) for tl in tlens]
idx = delta.index(min(delta))
best = subtrees[idx]
sys.stderr.write(f"Selected subtree with length {best.total_branch_length(): .3f} "
                 f"and {best.count_terminals()} tips\n")

# export sequences in this subtree
unique = {}
for tip in best.get_terminals():
    seq = fasta[tip.name]
    if seq in unique:
        # skip redundant sequences
        unique[seq].append(tip.name)
    else:
        unique.update({seq: [tip.name]})
        args.outfile.write(f">{tip.name}\n{seq}\n")

args.outfile.close()

if args.csvfile:
    for seq, labels in unique.items():
        if len(labels) == 1:
            continue
        key = labels[0]
        for label in labels[1:]:
            args.csvfile.write(f"{key},{label}\n")
    args.csvfile.close()
