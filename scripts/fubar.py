import argparse
import subprocess
import tempfile
import json
import csv
import os
import sys


description = """
Wrapper for the fast unconstrained Bayesian approximation 
(FUBAR) method in HyPhy. 
"""

def parse_args():
    """ Command line interface """
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("cds", type=str,
                        help="Path to codon alignment")
    parser.add_argument("prefix", type=str,
                        help="Path to write outputs as .fubar.csv and "
                        ".fubar.json")
    parser.add_argument("--hyphy", type=str, default="hyphy", 
                        help="Path to HyPhy executable")
    parser.add_argument("--ft2", type=str, default="fasttree", 
                        help="Path to fasttree executable")
    return parser.parse_args()


def clean_names(aln, hyphy_bin="hyphy", code="Universal", keep="No/No", stderr=None):
    """
    Clean stop codons and sequence names
    
    :param aln:  str, path to FASTA
    :param hyphy_bin:  str, path to HyPhy executable
    :param code:  str, genetic code (default: 'Universal')
    :param keep:  str, response to query "Keep all sequences and sites?"
    :param stderr:  for subprocess call, default to None
    :return:  str, path to temporary file with cleaned sequences
    """
    clean_file = tempfile.NamedTemporaryFile(delete=False)
    subprocess.run([
        hyphy_bin, 'cln', code, aln, keep, 
        clean_file.name  # destination file
        ], stdout=subprocess.PIPE, stderr=stderr)
    clean_file.close()
    return clean_file.name


def fasttree(infile, ft2_bin="fasttree", stderr=None):
    """
    Wrapper for fasttree

    :param infile:  str, path to input alignment file
    :param ft2_bin:  str, path to Fasttree binary
    :param stderr:  for subprocess call (default: None)
    :return:  str, path to temporary file with Newick tree string
    """
    tree_file = tempfile.NamedTemporaryFile(delete=False)
    subprocess.run([
        ft2_bin, '-nt',
        "-out", tree_file.name,
        infile
        ], stdout=subprocess.PIPE, stderr=stderr)
    tree_file.close()
    return tree_file.name


def fubar(aln, json_file, hyphy_bin="hyphy", ft2_bin="fasttree", verbose=False):
    """    
    :param aln: str, path to codon alignment
    :param json_file:  str, path to write JSON file
    :param hyphy_bin:  str, path to HYPHY executable
    :param ft2_bin:  str, path to Fasttree2 executable
    :param verbose:  bool, if True, stream messages to console
    :return: JSON output from FUBAR as a dict
    """
    stderr = None if verbose else subprocess.DEVNULL
    
    clean_file = clean_names(aln, hyphy=hyphy_bin, stderr=stderr)
    tree_file = fasttree(clean_file, ft2_bin=ft2_bin, stderr=stderr)
    
    # Run FUBAR
    result = subprocess.run([
        'hyphy', 'fubar', clean_file, tree_file,
        '--cache', '/dev/null', 
        '--output', json_file
        ], stdout=subprocess.PIPE, stderr=stderr)
    
    if verbose:
        for line in result.stdout:
            sys.stdout.write(line.decode('utf-8'))
    
    # load JSON into Python
    with open(json_file) as fp:
        results = json.load(fp)

    # clean up temp files
    os.remove(clean_file)
    os.remove(tree_file)
    
    return results
    

if __name__ == "__main__":
    args = parse_args()
    json_file = args.prefix + '.fubar.json'
    csv_file = args.prefix + '.fubar.csv'
    results = fubar(args.cds, json_file, hyphy_bin=args.hyphy, ft2_bin=args.ft2)
    
    # write MLE output as CSV
    with open(csv_file, 'w') as handle:
        writer = csv.writer(handle)
        writer.writerow(['pos'] + [label for label, _ in results["MLE"]['headers']])
        for i, row in enumerate(results["MLE"]['content']['0']):
            writer.writerow([i+1] + row[:-2])  # omit trailing "0,0"
