from Bio import SeqIO
import argparse
import sys
import os  

description = """
Generate a consensus sequence from a set of aligned sequences.
While this is possible in BioPython, it is rather convoluted.
"""

mixture_dict = {'W': 'AT', 'R': 'AG', 'K': 'GT', 'Y': 'CT', 'S': 'CG',
                'M': 'AC', 'V': 'AGC', 'H': 'ATC', 'D': 'ATG',
                'B': 'TGC', 'N': 'ATGC', '-': 'ATGC'}
ambig_dict = dict(("".join(sorted(v)), k) for k, v in mixture_dict.items())

def get_columns(aln, is_nuc=False):
    """ 
    Count frequencies of characters in each position 
    :param aln:  Bio.AlignIO object
    :param is_nuc:  bool, resolve nucleotide mixtures
    :return:  list of dicts (char, count) for all alignment positions
    """
    columns = None
    seqlen = None
    for record in aln:
        # handle first record
        if not columns:
            columns = [{nt: 1} for nt in record.seq]
            seqlen = len(record.seq)
            continue
        
        if len(record.seq) != seqlen:
            sys.stderr.write("ERROR: sequences are not the same length\n")
            sys.exit(1)
        
        for i, nt in enumerate(record.seq):
            if is_nuc and nt in mixture_dict:
                # count fractional states
                chars = mixture_dict[nt]
                for ch in chars:
                    if ch not in columns[i]:
                        columns[i].update({ch: 0})
                    columns[i][ch] += 1./len(chars)
            else:            
                if nt not in columns[i]:
                    columns[i].update({nt: 0})
                columns[i][nt] += 1
    return columns

def conseq(columns, thresh=0.5, is_nuc=False):
    """
    Generate consensus sequence from position-specific frequencies
    :param columns:  list of dicts from get_columns()
    :param thresh:  float, minimum frequency for character to be included in mixture
    :param is_nuc:  bool, if True then handle multiple states as a nucleotide mixture
    :return:  str
    """
    seq = ''
    for col in columns:
        if not col:  # Skip empty columns
            seq += '-'
            continue
        
        total = sum(col.values())
        if total == 0:  # Skip columns with no valid characters
            seq += '-'
            continue
        
        intermed = [(count/total, nt) for nt, count in col.items()]
        intermed.sort(reverse=True)
        states = [intermed[0][1]]  # always store most frequent character
        for freq, nt in intermed[1:]:
            if freq < thresh:
                break
            states.append(nt)  # mixtures
            
        if len(states) == 1:
            seq += states
        else:
            if not states:  # Handle empty states
                seq += '-'
            else:
                # multiple characters above threshold frequency
                if is_nuc:
                    key = ''.join(sorted(states))
                    seq += ambig_dict[key]  # convert nucleotide mixture
                else:
                    seq += states[0]  # append most frequent character
    return seq

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("infile", type=argparse.FileType('r'),
                        help="input, path to file containing aligned sequences")
    parser.add_argument("-f", "--format", default='fasta',
                        help="option, alignment file format (default: 'fasta')")
    parser.add_argument("-n", "--nuc", action="store_true",
                        help="Use nucleotide mixtures")
    parser.add_argument("-t", "--thresh", type=float, default=0.5,
                        help="If minority character is above this value (default 0.5), "
                        "then report '?' or mixture if '-n' is set.")
    args = parser.parse_args()
    
    try:
        aln = list(SeqIO.parse(args.infile, format=args.format))
        if not aln:
            sys.stderr.write("ERROR: Input file is empty or not in the correct format.\n")
            sys.exit(1)
        
        columns = get_columns(aln, is_nuc=args.nuc)
        if columns is None:
            sys.stderr.write("ERROR: No columns were generated. Check the input file.\n")
            sys.exit(1)
        
        cseq = conseq(columns, thresh=args.thresh, is_nuc=args.nuc)
        
        # get the names and extention 
        input_filename = os.path.splitext(os.path.basename(args.infile.name))[0]
        
        # write the consensus and save the file 
        sys.stdout.write(f">{input_filename}_consensus\n")  # Header
        sys.stdout.write(f"{cseq}\n")  # Secuencia consenso
    except Exception as e:
        sys.stderr.write(f"ERROR: {str(e)}\n")
        sys.exit(1)
