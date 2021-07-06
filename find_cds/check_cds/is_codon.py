"""
iterate through directory of fasta files 
summarizing how many of the fasta nuc seq entries are divisibile by 3
Summay of each file in a row on the csv outfile 
"""
import os

#directory = "/home/sareh/surfaces/find_cds/cut_cds/"
#out_path = "/home/sareh/surfaces/find_cds/check_cds/is_codon.csv"

directory = "/home/sareh/surfaces/find_cds/cut_cds"
out_path = "/home/sareh/surfaces/find_cds/check_cds/is_codon.csv"

def iter_fasta(handle):
    """
    Parse open file as FASTA.  Returns a generator
    of handle, sequence tuples.
    :param handle:  open stream to FASTA file in read mode
    :yield tuples, (header, sequence)
    """
    h, sequence = None, ''
    for line in handle:
        if line.startswith('>'):
            if len(sequence) > 0:
                yield h, sequence
                sequence = ''
            h = line.lstrip('>').rstrip()
        else:
            sequence += line.strip().upper()
    yield h, sequence


def main ():
    #args = parse_args()
    count = 0
    total_not_tripple=0
    file_count = 0 
    seq_count = 0
    #file_path = "/Users/sareh/Desktop/cutter_NC_001436_NP_057863.1"
    with open (out_path, 'w') as outfile:
        outfile.write("file_name,count_not_tripple,lst_not_tripple")
        for file in os.listdir(directory):
            not_tripple = 0
            file_count += 1 
            lst_not_tripple = []
            with open(file, 'r') as f:
                for h,sequence in iter_fasta(f):
                    seq_count += 1 
                    divisible = len(sequence)% 3
                    if divisible != 0:
                        lst_not_tripple.append(h)
                        #print("{},{}".format(divisible,h))
                        not_tripple += 1
                outfile.write("{},{},{}".format(file,not_tripple,lst_not_tripple))

    print(total_not_tripple)
    print(file_count) 
    print(seq_count)

if __name__ == '__main__':
    main()
