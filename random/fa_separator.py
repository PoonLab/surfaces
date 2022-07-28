
#!/usr/bin python3

import os
import argparse

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--ref', type=str,
                       help='fasta file')
    parser.add_argument('--dir', type=str,
                       help='dir to write out every file')
    return parser.parse_args()


def convert_fasta (handle):
    """
    takes in handel to an open fasta file
    """
    result = []
    sequence = ''
    # handle = open(file,'r')
    for line in handle:
        if line.startswith('$'): # skip header line
            continue
        elif line.startswith('>') or line.startswith('#'):
            if len(sequence) > 0:
                result.append([h,sequence])
                sequence = ''   # reset
            h = line.strip('>#\n')
        else:
            sequence += line.strip('\n').upper()

    result.append([h,sequence]) # handle last entry
    return result

def main():

    args = parse_args()

    out_dir = args.dir
    ref_file = args.ref
    ref_handle = open(ref_file, "r")

    ref_tuple = convert_fasta(ref_handle)
    for ref_h, ref_seq in ref_tuple:
        ref_accn = ref_h.split(",")[0]
        out_path = os.path.join(out_dir,ref_accn)
        out_handle = open(out_path, "w")
        out_handle.write(">{}\n{}\n".format(ref_h, ref_seq))

if __name__ == "__main__":
    main()


