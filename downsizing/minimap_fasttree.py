import os
import re
import argparse
import subprocess

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
    return parser.parse_args()

def main():
    args = parse_args()
    count = 0

    for file in os.listdir(args.indir):
        in_path = os.path.join(args.indir, file)
        out_path = os.path.join(args.outdir, file)
        count += 1
        if count % total_number != my_number:
            continue

        out_handle = open(out_path, "w")

        #2. fasttree -nt -gtr -gamma minimap_alignment.fasta > fasttree.tre
        #os.system('fasttree' '-nt' '-gtr' '-gamma' in_path '>' out_path)
        #subprocess.call(cmd_fasttree, stdout=out_handle, shell=True)
        #cmd_fasttree = ['fasttree', '-nt', '-gtr', '-gamma', str(in_path), '>', str(out_path)]

        cmd_fasttree = ['fasttree', '-nt', '-gtr', '-gamma', str(in_path)]
        ftout = subprocess.run(cmd_fasttree, stdout=subprocess.PIPE)

        #writing it out
        a = str(ftout.stdout)
       
        b = a[2:] #remove: b'
        
        c = b[:-3] #remove: '\n
        
        out_handle.write(c)

if __name__ == "__main__":
    main()
