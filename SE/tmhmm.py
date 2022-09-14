import os
import biolib
import argparse
from IPython.display import Image

#out_dir = "/home/sareh/all/tmhmm"

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--indir', type=str,
                       help='directory of aa fasta files')
    parser.add_argument('--outdir', type=str,
                        help='tmhmm output')
    return parser.parse_args()

args = parse_args()

for accn in os.listdir(args.indir):
    infile = os.path.join(args.indir, accn)
    print(infile)
    outfile = os.path.join(args.outdir, accn)
    print(outfile)

    deeptmhmm = biolib.load('DTU/DeepTMHMM')
    # Run DeepTMHMM
    tmhmm = "--fasta " + infile
    print(tmhmm)
    deeptmhmm_res = deeptmhmm.cli(args=tmhmm)

    print(deeptmhmm_res)

    # Save the results
    deeptmhmm_res.save_files(outfile)

# Error with the image file
#png_out = os.path.join(out_dir,"plot",accn)
#Image(filename=png_out)
