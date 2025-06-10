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
args = parser.parse_args()

# import data from FASTA file
records = SeqIO.parse(args.fasta, "fasta")
fasta = {}
for record in records:
    fasta.update({record.description: record.seq})

# import tree and check labels against FASTA
tree = Phylo.read(args.tree, 'newick')
tips = [tip.name for tip in tree.get_terminals()]
for tip in tips:
    if tip not in fasta:
        sys.stderr(f"ERROR: Could not locate {tip} in FASTA!\n")
        sys.exit(1)
if args.branch:
    assert args.branch in tips, f"--branch {args.branch} is not in tree!"

# by default, export only sequences that appear in the tree
if args.target is None or args.branch is None:
    sys.stderr.write("No target set by user, writing all tips in tree to FASTA by default...\n")
    for tip in tips:
        args.outfile.write(f">{tip}\n{fasta[tip]}\n")
    args.outfile.close()
    sys.exit()


assert args.target > 0, "--target must be greater than zero!"
assert args.target < tree.total_branch_length(), "--target exceeds length of entire tree!"

# otherwise traverse from user-specified tip toward root
nodes = tree.get_path(args.branch)
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
for tip in best.get_terminals():
    args.outfile.write(f">{tip.name}\n{fasta[tip.name]}\n")

args.outfile.close()
