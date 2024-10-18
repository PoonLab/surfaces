from prunetree import prunetree
import argparse
from glob import glob
import sys
import os

parser = argparse.ArgumentParser("Batch processing with prunetree")
parser.add_argument("csvfile", type=argparse.FileType('r'),
                    help="CSV output from step4_filter.R")
parser.add_argument("glob", type=str, 
                    help="Quoted glob for Newick files, e.g., \"data/Virus_*_step3.nwk\"")
parser.add_argument("--overwrite", action="store_true", 
                    help="Option to overwrite output files")
args = parser.parse_args()

files = glob(args.glob)
cutoffs = {}
for line in args.csvfile:
    prot, ntips = line.strip().split(',')
    cutoffs.update({prot: int(ntips)})

for f in files:
    if not f.endswith("_step3.nwk"):
        sys.stderr.write("ERROR: filename must end with _step3.nwk suffix\n")
        sys.stderr.write(f"Failed {f}\n")
        sys.exit()
    prot = os.path.basename(f).split('_')[1]
    if prot not in cutoffs:
        sys.stderr.write(f"ERROR: unrecognized protein {prot} in filename {f}\n")
    
    # prepare output filenames
    seqfile = f.replace(".nwk", ".fasta")
    if not os.path.exists(seqfile):
        sys.stderr.write(f"ERROR: Failed to locate seqfile {seqfile}\n")
        sys.exit()
    
    outfile = f.replace("_step3.nwk", "_step4.fasta")
    if os.path.exists(outfile) and not args.overwrite:
        sys.stderr.write(f"ERROR: output file {outfile} exists (use --overwrite?)\n")
        sys.exit()
    outfile = open(outfile, 'w')
    
    csvfile = f.replace("_step3.nwk", "_step4.labels.csv")
    if os.path.exists(csvfile) and not args.overwrite:
        sys.stderr.write(f"ERROR: output file {csvfile} exists (use --overwrite?)\n")
        sys.exit()
    csvfile = open(csvfile, 'w')
    
    prunetree(f, mode='ntips', target=cutoffs[prot], seq=seqfile, 
              outfile=outfile, csvfile=csvfile)
