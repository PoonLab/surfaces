from fubar import fubar
from glob import glob
from mpi4py import MPI
import sys
import os

nprocs = MPI.COMM_WORLD.Get_size()
my_rank = MPI.COMM_WORLD.Get_rank()

if nprocs == 1:
    sys.stderr.write("Detected only one process - did you use mpirun? Exiting.\n")
    sys.exit()


#files = glob("data/5_alt/*_step5.fasta")
files = glob("data/8_fubar/*_step7.fasta")

for i, f in enumerate(files):
    if i % nprocs != my_rank:
        continue
    sys.stderr.write(f"{my_rank}/{nprocs} starting job {i} of {len(files)}: {os.path.basename(f)}\n")
    #out = f.replace("5_alt", "6_fubar").replace(".fasta", ".fubar.json")
    out = f.replace(".fasta", ".fubar.json")
    fubar(aln=f, json_file=out, verbose=False)

