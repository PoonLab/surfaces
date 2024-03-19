"""
To estimate ideal tree length, use indelible to simulate sequences under different scaling factors
Use this simulated sequences to measure selection and estimate differences between simualated and calculated omegas
"""

import subprocess
import argparse
import sys
from pathlib import Path
from selection_by_cluster import run_selection_pipeline


def get_args(parser):
# Arguments of input and output
    parser.add_argument('tree', help='Input the tree')
    parser.add_argument('-c', '--control', help='Output file path', 
                        default = 'control.txt')
    parser.add_argument('-p', '--path', help='Path to store simulation information', 
                        default = 'indelible.out')
    parser.add_argument('-r', '--replicates', default = 1, 
                        help='Run indelible from control file')
    return parser.parse_args()


tree_lengths = [0.1, 0.25, 0.5, 0.75, 1, 1.5, 2, 5, 10, 15, 20]

if __name__=="__main__":
    parser = argparse.ArgumentParser(
        description="Create control file to run indelible"
    )

    args = get_args(parser)

    for i in tree_lengths:
        sim_name = Path(args.path, f'sim_{i}')
        # Create control files and run indelible
        print(sim_name)

        subprocess.run(['python3', 'create_control_indelible.py', args.tree, 
                        '-o', str(sim_name), 
                        '-rep', str(args.replicates),
                        '-s', str(i), '-r'], check=True)
        
        alignment = f'{sim_name}.fas'
        cleaned_file = f'{sim_name}.cleaned.fa'
        tree_file = f'{sim_name}.cleaned.nwk'
        out_name = f'{sim_name}.FUBAR.json'

        # Measure FUBAR
        run_selection_pipeline(alignment, cleaned_file, tree_file, out_name)
        # sys.exit()




