import argparse
import os
import subprocess
import sys
from glob import glob


description = """
Run a HyPhy selection method on simulated tree-alignment pairs that have
already been generated and cleaned.
"""


def parse_args():
    """Command line interface."""
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        "indir", type=str,
        help="Path to directory containing simulation files"
    )
    parser.add_argument(
        "--method", choices=["slac", "fel", "fubar"], default="slac",
        help="HyPhy method to run"
    )
    parser.add_argument(
        "--outdir", type=str, default=None,
        help="Name of output subdirectory to create inside indir"
    )
    parser.add_argument(
        "--hyphy", type=str, default="hyphy",
        help="Path to HyPhy executable"
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


def run_method(aln, tree, json_file, method="slac", hyphy_bin="hyphy",
               verbose=False):
    """
    Run a HyPhy selection method on a cleaned alignment and matching tree.
    """
    stderr = None if verbose else subprocess.DEVNULL
    cmd = [
        hyphy_bin, method,
        "--alignment", aln,
        "--tree", tree,
        "--output", json_file
    ]

    if method in {"slac", "fel"}:
        cmd.extend(["--branches", "All"])

    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=stderr)

    if verbose:
        sys.stdout.write(result.stdout.decode("utf-8"))

    return result.returncode


if __name__ == "__main__":
    args = parse_args()

    indir = os.path.abspath(args.indir)
    method_label = args.method.upper()
    outdir_name = args.outdir if args.outdir else f"{args.method}_outputs"
    outdir = os.path.join(indir, outdir_name)
    os.makedirs(outdir, exist_ok=True)

    alignments = sorted(glob(os.path.join(indir, "*.cleaned.fa")))

    n_success = 0
    n_existing = 0
    n_incomplete = 0
    n_failed = 0

    for aln in alignments:
        prefix = os.path.basename(aln).replace(".cleaned.fa", "")
        tree = os.path.join(indir, prefix + ".cleaned.nwk")
        outfile = os.path.join(outdir, prefix + f".{method_label}.json")

        if not os.path.exists(tree):
            sys.stderr.write(f"...skipping incomplete pair for {prefix}\n")
            n_incomplete += 1
            continue

        if os.path.exists(outfile) and not args.overwrite:
            sys.stderr.write(
                f"...skipping existing output {os.path.basename(outfile)}\n"
            )
            n_existing += 1
            continue

        sys.stderr.write(f"...running {method_label} on {prefix}\n")
        sys.stderr.flush()

        code = run_method(
            aln=aln,
            tree=tree,
            json_file=outfile,
            method=args.method,
            hyphy_bin=args.hyphy,
            verbose=args.verbose
        )

        if code == 0 and os.path.exists(outfile):
            n_success += 1
        else:
            sys.stderr.write(f"...{method_label} failed for {prefix}\n")
            n_failed += 1

    sys.stderr.write(
        f"Completed: {len(alignments)} alignments scanned, "
        f"{n_success} succeeded, "
        f"{n_existing} skipped existing, "
        f"{n_incomplete} incomplete, "
        f"{n_failed} failed\n"
    )
