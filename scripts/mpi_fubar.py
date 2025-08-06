from fubar import fubar
from glob import glob
from mpi4py import MPI
import sys
import os
import argparse

desc = """
Wrapper to run the FUBAR method in HyPhy within an MPI parallel
computing environment.
"""

nprocs = MPI.COMM_WORLD.Get_size()
my_rank = MPI.COMM_WORLD.Get_rank()

parser = argparse.ArgumentParser(description=desc)
parser.add_argument("glob", type=str, help="UNIX glob for FASTA files to process. "
                    "This argument MUST be enclosed in double quotes, e.g., "
                    "\"data/*_step5.fasta\"")
parser.add_argument("--dry-run", action="store_true", 
                    help="List files for processing without actually running FUBAR.")
parser.add_argument("--overwrite", action="store_true",
                    help="Replace any existing output files.")
args = parser.parse_args()

if nprocs == 1:
    sys.stderr.write("Detected only one process - did you use mpirun? Exiting.\n")
    sys.exit()


files = glob(args.glob)

for i, f in enumerate(files):
    if i % nprocs != my_rank:
        continue
    
    out = f.replace(".fasta", ".fubar.json")
    if os.path.exists(out):
        if args.dry_run:
            if args.overwrite:
                sys.stderr.write(f"...would OVERWRITE {out}\n")
            else:
                sys.stderr.write(f"...would skip existing file {out}\n")
            continue
        else:
            # actual run!
            if not args.overwrite:
                sys.stderr.write(f"...skipping existing file {out}\n")
                continue  # go to next file
    
    if args.dry_run:
        sys.stderr.write(f"...would write to {out}\n")
        continue
    
    sys.stderr.write(f"{my_rank}/{nprocs} starting job {i} of {len(files)}: "
                     f"{os.path.basename(f)}\n")
    sys.stderr.flush()
    fubar(aln=f, json_file=out, verbose=False)

