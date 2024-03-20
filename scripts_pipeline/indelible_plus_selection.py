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
    parser.add_argument('-r', '--replicates', default = 1, type=int,
                        help='Run indelible from control file')
    return parser.parse_args()


# tree_lengths = [ 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09,
#                  0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
#                  1.0, 2.0, 3.0, 4.0, 5.0, 10.0, 15.0, 20.0, 30.0]

tree_lengths = [ 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09]

if __name__=="__main__":
    parser = argparse.ArgumentParser(
        description="Create control file to run indelible"
    )

    args = get_args(parser)

    for i in tree_lengths:
        
        for j in list(range(args.replicates)):
            sim_name = Path(args.path, f'sim_{i}_{j}_')
            
            print(f'Replicate: {j}; Tree lenght: {i}')
            #Create control files and run indelible

            subprocess.run(['python3', 'create_control_indelible.py', args.tree, 
                            '-o', str(sim_name), 
                            # '-rep', str(args.replicates),  # creates complicated outputs
                            '-s', str(i), '-r'], check=True)
            
            alignment = f'{sim_name}.fas'
            cleaned_file = f'{sim_name}.cleaned.fa'
            tree_file = f'{sim_name}.cleaned.nwk'
            out_name = f'{sim_name}.FUBAR.json'

            # Measure FUBAR
            run_selection_pipeline(alignment, cleaned_file, tree_file, out_name)
            # sys.exit()




