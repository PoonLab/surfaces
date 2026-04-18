import argparse
import subprocess
import sys
import os
from glob import glob


description = """
Run HyPhy SLAC on simulated tree-alignment pairs that have already been
generated and cleaned.
"""


def parse_args():
    """Command line interface."""
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        "indir", type=str,
        help="Path to directory containing simulation files"
    )
    parser.add_argument(
        "--outdir", type=str, default="slac_outputs",
        help="Name of output subdirectory to create inside indir"
    )
    parser.add_argument(
        "--hyphy", type=str, default="hyphy",
        help="Path to HyPhy executable"
    )
    parser.add_argument(
        "--code", type=str, default="Universal",
        help="Genetic code to pass to HyPhy SLAC"
    )
    parser.add_argument(
        "--branches", type=str, default="All",
        help="Branches setting for HyPhy SLAC"
    )
    parser.add_argument(
        "--samples", type=int, default=0,
        help="Number of samples for ancestral reconstruction uncertainty"
    )
    parser.add_argument(
        "--pvalue", type=float, default=0.1,
        help="p-value reporting threshold"
    )
    parser.add_argument(
        "--overwrite", action="store_true",
        help="Overwrite existing output JSON files"
    )
    parser.add_argument(
        "--verbose", action="store_true",
        help="Display HyPhy messages"
    )
    return parser.parse_args()


def run_slac(aln, tree, json_file, hyphy_bin="hyphy", code="Universal",
             branches="All", samples=0, pvalue=0.1, verbose=False):
    """
    Run HyPhy SLAC on a cleaned alignment and matching tree.
    """
    stderr = None if verbose else subprocess.DEVNULL
    result = subprocess.run([
        hyphy_bin, "slac",
        "--code", code,
        "--alignment", aln,
        "--tree", tree,
        "--branches", str(branches),
        "--samples", str(samples),
        "--pvalue", str(pvalue),
        "--output", json_file
    ], stdout=subprocess.PIPE, stderr=stderr)

    if verbose:
        sys.stdout.write(result.stdout.decode("utf-8"))

    return result.returncode


if __name__ == "__main__":
    args = parse_args()

    indir = os.path.abspath(args.indir)
    outdir = os.path.join(indir, args.outdir)
    os.makedirs(outdir, exist_ok=True)

    alignments = sorted(glob(os.path.join(indir, "*.cleaned.fa")))

    n_success = 0
    n_existing = 0
    n_incomplete = 0
    n_failed = 0

    for aln in alignments:
        prefix = os.path.basename(aln).replace(".cleaned.fa", "")
        tree = os.path.join(indir, prefix + ".cleaned.nwk")
        outfile = os.path.join(outdir, prefix + ".SLAC.json")

        if not os.path.exists(tree):
            sys.stderr.write(f"...skipping incomplete pair for {prefix}\n")
            n_incomplete += 1
            continue

        if os.path.exists(outfile) and not args.overwrite:
            sys.stderr.write(f"...skipping existing output {os.path.basename(outfile)}\n")
            n_existing += 1
            continue

        sys.stderr.write(f"...running SLAC on {prefix}\n")
        sys.stderr.flush()

        code = run_slac(
            aln=aln,
            tree=tree,
            json_file=outfile,
            hyphy_bin=args.hyphy,
            code=args.code,
            branches=args.branches,
            samples=args.samples,
            pvalue=args.pvalue,
            verbose=args.verbose
        )

        if code == 0 and os.path.exists(outfile):
            n_success += 1
        else:
            sys.stderr.write(f"...SLAC failed for {prefix}\n")
            n_failed += 1

    sys.stderr.write(
        f"Completed: {len(alignments)} alignments scanned, "
        f"{n_success} succeeded, "
        f"{n_existing} skipped existing, "
        f"{n_incomplete} incomplete, "
        f"{n_failed} failed\n"
    )
