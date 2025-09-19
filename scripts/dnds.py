import argparse
import json
import subprocess
import sys
import os
import tempfile
from glob import glob
import csv

from fubar import clean_names, fasttree


description = """
Wrapper for the fast unconstrained Bayesian approximation 
(FUBAR) method in HyPhy. 
"""

def slac(aln, code='Universal', branches='All', samples=0, pval=0.1, 
         hyphy_bin="hyphy", ft2_bin="fasttree", verbose=False):
    """    
    Wrapper for HyPhy SLAC method

    :param aln:  str, path to codon alignment
    :param code:  str, genetic code (default: 'Universal')
    :param branches:  str, branches to test (default: 'All')
    :param samples:  int, number of samples to measure uncertainty in 
                     ancestral reconstruction (default: 0)
    :param pval:  float, p-value threshold for reporting (default: 0.1)
    :param hyphy_bin:  str, path to HYPHY executable
    :param ft2_bin:  str, path to Fasttree2 executable
    :param verbose:  bool, if True, stream messages to console
    :return: JSON output from FUBAR as a dict
    """
    stderr = None if verbose else subprocess.DEVNULL
    
    clean_file = clean_names(aln, hyphy_bin=hyphy_bin, stderr=stderr)
    tree_file = fasttree(clean_file, ft2_bin=ft2_bin, stderr=stderr)
    json_file = tempfile.NamedTemporaryFile(delete=False)

    # Run FUBAR
    result = subprocess.run([
        'hyphy', 'slac', 
        '--code', code,
        '--alignment', clean_file, 
        '--tree', tree_file,
        '--branches', str(branches),
        '--samples', str(samples),
        '--pvalue', str(pval),
        '--output', json_file.name
        ], stdout=subprocess.PIPE, stderr=stderr)
    
    if verbose:
        for line in result.stdout.decode('utf-8'):
            try:
                sys.stdout.write(line)
            except AttributeError:
                print(result.stdout)
                raise
    
    # load JSON into Python
    with open(json_file.name) as fp:
        results = json.load(json_file)

    # clean up temp files
    os.remove(clean_file)
    os.remove(tree_file)
    os.remove(json_file.name)
    
    return results
    

def mean_dnds(aln, code='Universal', branches='All', samples=0, pval=0.1, 
         hyphy_bin="hyphy", ft2_bin="fasttree", verbose=False, ambig='AVERAGED'):
    """ 
    Wrapper to extract mean dN/dS ratio from SLAC analysis 
    
    :param aln:  str, path to codon alignment
    :param code:  str, genetic code (default: 'Universal')
    :param branches:  str, branches to test (default: 'All')
    :param samples:  int, number of samples to measure uncertainty in 
                     ancestral reconstruction (default: 0)
    :param pval:  float, p-value threshold for reporting (default: 0.1)
    :param hyphy_bin:  str, path to HYPHY executable
    :param ft2_bin:  str, path to Fasttree2 executable
    :param verbose:  bool, if True, stream messages to console
    :param ambig:  str, AVERAGED (default) or RESOLVED.  If AVERAGED, then 
                   all possible resolutions of an ambiguous character contribute 
                   to calculating dN and dS.
    :return:  tuple, (global dN-dS, number of sites under negative (purifying) selection,
              number of sites under positive (diversifying) selection, and total number 
              of codon sites)
    """
    res = slac(aln, code=code, branches=branches, samples=samples, pval=pval,
               hyphy_bin=hyphy_bin, ft2_bin=ft2_bin, verbose=verbose)
    global_dnds = res['fits']['Global MG94xREV']['Rate Distributions'][
        'non-synonymous/synonymous rate ratio for *test*'][0][0]
    
    n_pos, n_neg, n_cod = 0, 0, 0
    for row in res['MLE']['content']['0']['by-site'][ambig]:
        es, en, s, n, ps, ds, dn, dnds, p_pos, p_neg, tbl = row
        n_cod += 1
        if p_pos < pval:
            n_pos += 1
        if p_neg < pval:
            n_neg += 1
    
    return global_dnds, n_pos, n_neg, n_cod


if __name__ == "__main__":
    # Command line interface
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("glob", type=str,
                        help="UNIX glob for input files")

    parser.add_argument("--hyphy", type=str, default="hyphy", 
                        help="Path to HyPhy executable")
    parser.add_argument("--ft2", type=str, default="fasttree", 
                        help="Path to fasttree executable")
    parser.add_argument("--pval", type=float, default=0.1,
                        help="p-value reporting threshold (default: 0.1)")
    parser.add_argument("--verbose", action="store_true", 
                        help="option, display messages")
    parser.add_argument("-o", "--outfile", type=argparse.FileType('w'),
                        default=sys.stdout, 
                        help="File to write CSV output (default: stdout)")
    args = parser.parse_args()

    files = glob(args.glob)
    writer = csv.writer(args.outfile)
    writer.writerow(['virus', 'protein', 'dnds', 'n.pos', 'n.neg', 'codons'])

    for path in files:
        fn = os.path.basename(path)
        tokens = fn.split('.')[0].split('_')
        virus = tokens[0]
        protein = ' '.join(tokens[1:-1])

        dnds, n_pos, n_neg, total = mean_dnds(path, hyphy_bin=args.hyphy, ft2_bin=args.ft2, 
                                    verbose=args.verbose, pval=args.pval)

        writer.writerow([virus, protein, dnds, n_pos, n_neg, total])
        args.outfile.flush()
 

