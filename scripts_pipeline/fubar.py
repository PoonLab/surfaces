import argparse
import subprocess
import tempfile
import json
import csv
import os


description = """
Wrapper for the fast unconstrained Bayesian approximation 
(FUBAR) method in HyPhy. 
"""

def parse_args():
    """ Command line interface """
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("cds", type=str,
                        help="Path to codon alignment")
    parser.add_argument("csvfile", type=argparse.FileType('w'),
                        help="Path to write FUBAR output as CSV")
    parser.add_argument("--hyphy", type=str, default="hyphy", 
                        help="Path to HyPhy executable")
    parser.add_argument("--ft2", type=str, default="fasttree", 
                        help="Path to fasttree executable")
    return parser.parse_args()


def fubar(aln, hyphy_bin="hyphy", ft2_bin="fasttree"):
    """    
    :param aln: str, path to codon alignment
    :return: JSON output from FUBAR as a dict
    """
    
    # Clean stop codons and sequence names
    clean_file = tempfile.NamedTemporaryFile(delete=False)
    subprocess.run([
        hyphy_bin, 'cln',
        'Universal',  # genetic code 
        aln, 'No/No',  # Keep all sequences and sites
        clean_file.name  # destination file
        ], stdout=subprocess.PIPE)
    clean_file.close()

    # Regenerate tree, sequence names may have changed
    tree_file = tempfile.NamedTemporaryFile(delete=False)
    subprocess.run([
        ft2_bin, '-nt',
        "-out", tree_file.name,
        clean_file.name
        ], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    tree_file.close()
    
    # Run FUBAR
    json_file = tempfile.NamedTemporaryFile(delete=False)
    p = subprocess.run([
        'hyphy', 'fubar', clean_file.name, tree_file.name,
        '--cache', '/dev/null', 
        '--output', json_file.name
        ], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    json_file.close()
    
    # load JSON into Python
    with open(json_file.name) as fp:
        results = json.load(fp)

    # clean up temp files
    os.remove(clean_file.name)
    os.remove(tree_file.name)
    os.remove(json_file.name)
    
    return results
    

if __name__ == "__main__":
    args = parse_args()
    results = fubar(args.cds, hyphy_bin=args.hyphy, ft2_bin=args.ft2)
    
    # write MLE output as CSV
    writer = csv.writer(args.csvfile)
    writer.writerow(['pos'] + [label for label, _ in results["MLE"]['headers']])
    for i, row in enumerate(results["MLE"]['content']['0']):
        writer.writerow([i+1] + row[:-2])  # omit trailing "0,0"
    args.csvfile.close()
