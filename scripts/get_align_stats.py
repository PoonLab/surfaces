import subprocess
import tempfile
import os
from Bio import AlignIO, Phylo
from glob import glob
import argparse


def get_stats(fasta, hyphy_bin="hyphy", ft2_bin="fasttree", verbose=False):
    """
    :param fasta:  path to codon alignment (FASTA format)
    :param hyphy_bin:  path to HyPhy executable
    :param ft2_bin:  path to fasttree2 executable
    :param verbose:  bool, display messages
    """
    stderr = None if verbose else subprocess.DEVNULL
    
    # Clean stop codons and sequence names
    clean_file = tempfile.NamedTemporaryFile(delete=False)
    subprocess.run([
        hyphy_bin, 'cln',
        'Universal',  # genetic code 
        fasta, 'No/No',  # Keep all sequences and sites
        clean_file.name  # destination file
        ], stdout=subprocess.PIPE, stderr=stderr)
    clean_file.close()
    
    aln = AlignIO.read(clean_file, format='fasta')
    ncod = aln.get_alignment_length()//3
    nseq = len(aln)

    # Regenerate tree, sequence names may have changed
    tree_file = tempfile.NamedTemporaryFile(delete=False)
    subprocess.run([
        ft2_bin, '-nt',
        "-out", tree_file.name,
        clean_file.name
        ], stdout=subprocess.PIPE, stderr=stderr)
    tree_file.close()
    
    phy = Phylo.read(tree_file, 'newick')
    treelen = phy.get_total_branch_length()
    
    # clean up temp files
    os.remove(clean_file.name)
    os.remove(tree_file.name)

    return {'ncod': ncod, 'nseq': nseq, 'treelen': treelen}


parser = argparse.ArgumentParser()
files = glob()