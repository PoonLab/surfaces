"""
Create control file and run indelible 
"""

#import libraries 
import argparse
import subprocess
import re
import os


def get_args(parser):
# Arguments of input and output
    parser.add_argument('tree', help='Input the tree')
    parser.add_argument('-c', '--control', help='Output file path', default = 'control.txt')
    parser.add_argument('-o', '--out', help='Output name', default = 'indelible.out')
    parser.add_argument('-s', '--scaling_factor', type = int, default=20,
                            help='Scaling factor to define tree length')
    parser.add_argument('-r', '--run', action = 'store_true', help='Run indelible from control file')
    return parser.parse_args()

def create_control(control, prop, kappa, omegas, scaling_factor, tree_string, n_sites, outname):
    """
    Create control file to run indelible
    :param output: str, name of output file
    :param prop: list, proportions of omegas
    :param omegas: list, omega values 
    :param scaling_factor: int, tree length 
    :param tree_string: str, phylogenetic tree as read from newick file
    :param n_sites: int, number of sites in the sequence
    """
    # write minimal contents of INDELible control file
    with open(control, 'w') as handle:
        handle.write(f'[TYPE] CODON 1\n')
        handle.write(f'[SETTINGS]\n[output] FASTA\n[printrates] TRUE\n') 
        handle.write(f'[MODEL] M3\n[submodel] {kappa}\n')

        prop_string = ' '.join(map(lambda p: '%1.9f' % p, prop))
        omega_string = ' '.join(map(str, omegas))

        handle.write(f' {prop_string}\n {omega_string}\n')

        #handle.write('[indelmodel] NB 0.2 4\n')
        handle.write(f'[indelmodel] LAV 1.5 4\n')
        handle.write(f'[indelrate] 0.0\n')

        handle.write('[TREE] bigtree %s;\n' % tree_string.rstrip('0123456789.;:\n'))
        handle.write(f'[treelength] {scaling_factor}\n')
        handle.write('[PARTITIONS] partitionname\n')
        handle.write(f'  [bigtree M3 {n_sites}]\n')  # treename name rootlength
        #handle.write("  [bigtree M3 HIV-pol-KC169753.txt]\n")
        handle.write(f'[EVOLVE] partitionname 1 {outname} \n' ) 
        print(f'Control file at: {control}')
        print(f'output name: {outname}')

    handle.close()

def run_indelible(control_file):
    """Run indelible from control file"""
    
    subprocess.run(['indelible', control_file])

if __name__=="__main__":
    parser = argparse.ArgumentParser(
        description="Create control file to run indelible"
    )

    args = get_args(parser)
    
    # Define inputs
    outname = args.out
    control = args.control
    handle = open(args.tree, 'r')
    tree_string = handle.readline()
    handle.close()

    #print(tree_string)

    #filename = os.path.basename(args.file)
    #prefix = filename.split('.')[0]
    #print(prefix)

    # factor to stretch entire tree by (mutation rate)
    #global_scaling_factors = [1.0]
    #target_depth = 0.1  # maximum expected number of substitutions from tip to root
    #tree = Phylo.read(treefile, 'newick')
    #max_depth = max(tree.depths().values())
    #scaling_factor = target_depth / max_depth
    scaling_factor = args.scaling_factor

    #omegas = [0.01, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0]

    # Gamma(shape=1.5, rate=3), histogram with 50 breaks
    omegas = [
    0.0510204, 0.1530612, 0.2551020, 0.3571428, 0.4591836, 0.5612244, 0.6632652, 0.7653060, 0.8673468,
    0.9693876, 1.0714284, 1.1734692, 1.2755100, 1.3775508, 1.4795916, 1.5816324, 1.6836732, 1.7857140,
    1.8877548, 1.9897956, 2.0918364, 2.1938772, 2.2959180, 2.3979588, 2.4999996, 2.6020404, 2.7040812,
    2.8061220, 2.9081628, 3.0102036, 3.1122444, 3.2142852, 3.3163260, 3.4183668, 3.5204076, 3.6224484,
    3.7244892, 3.8265300, 3.9285708, 4.0306116, 4.1326524, 4.2346932, 4.3367340, 4.4387748, 4.5408156,
    4.6428564, 4.7448972, 4.8469380, 4.9489788, 5.0510196]

    prop = [
        1.063764e-01, 1.464867e-01, 1.401630e-01, 1.223917e-01, 1.023007e-01, 8.333079e-02, 6.673162e-02, 
        5.279576e-02, 4.139382e-02, 3.222714e-02, 2.495008e-02, 1.922789e-02, 1.476163e-02, 1.129629e-02,
        8.620596e-03, 6.562952e-03, 4.985996e-03, 3.780965e-03, 2.862470e-03, 2.163923e-03, 1.633688e-03,
        1.231906e-03, 9.279287e-04, 6.982665e-04, 5.249686e-04, 3.943506e-04, 2.960034e-04, 2.220245e-04,
        1.664245e-04, 1.246710e-04, 9.333909e-05, 6.984377e-05, 5.223631e-05, 3.904912e-05, 2.917804e-05,
        2.179306e-05, 1.627075e-05, 1.214322e-05, 9.059520e-06, 6.756626e-06, 5.037501e-06, 3.754636e-06,
        2.797655e-06, 2.084012e-06, 1.551998e-06, 1.155507e-06, 8.600995e-07, 6.400653e-07, 4.762155e-07
    ]

    n_sites = 1000 # codons
    kappa = 8.0 

    # Create file
    create_control(control, prop, kappa, omegas, scaling_factor, tree_string, n_sites, outname)
    if args.run:
        run_indelible(control)
 

