import os
import re
import argparse
import subprocess

mm2bin = "/home/sareh/2020/scripts/downsizing/minimap2.py"

# mpi
from mpi4py import MPI
the_world = MPI.COMM_WORLD
my_number = the_world.Get_rank()
total_number = the_world.Get_size()

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--indir', type=str,
                       help='IN: directory of input files')
    parser.add_argument('--outdir', type=str,
                       help='OUT: directory to output file')
    parser.add_argument('--ref', type=str,
                       help='REF: directory to refrence genome fasta file')
    return parser.parse_args()

def main():

    nthread = 5
    args = parse_args()
    count = 0

    for file in os.listdir(args.indir):
        in_path = os.path.join(args.indir, file)
        out_path = os.path.join(args.outdir, file)
        ref_path = os.path.join(args.ref, file)
        count += 1
        if count % total_number != my_number:
            #sleep(1)
            continue

        #arg = ("python3 minimap.py -o {} -ref {} -f -a".format(out_path, in_path,))
        #os.system(arg)

        #1. python3 minimap2.py -o minimap_alignment.fasta --ref ref.hiv.test.fa NC_001802.txt -f -a
        cmd_minimap = ['python3', mm2bin, '-o',str(out_path), str(in_path), '--ref', str(ref_path), '-f', '-a', '-t$
        subprocess.run(cmd_minimap)
        print(file)

if __name__ == "__main__":
    main()
