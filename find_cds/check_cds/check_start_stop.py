import re
import os


directory="/home/sareh/surfaces/find_cds/cut_cds"

start = ["ATG"] #
stop = ["TAA","TAG","TGA"]

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


def main():
    total = 0
    start_correct = 0
    start_wrong = 0
    lst_start=[]
    lst_stop=[]
    stop = ["TAA", "TAG", "TGA"]
    for f in (os.listdir(directory)):
        path = os.path.join(directory,f)
        with open(path, 'r') as file:
            for h, sequence in iter_fasta(file):
                total += 1
                #CHECKed RNA seq it's all T
                start_position=re.search("ATG",sequence)
                stop_Position1=re.search("TAA",sequence)
                stop_Position2=re.search ("TGA",sequence)
                stop_Position3=re.search ("TAG",sequence)
                # stop_patter = '(?<=ATG)(?:[ACGT]{3})*(TAA|TGA|TAG)'
                seq_len=len(sequence)
                lst_stop.append(sequence[seq_len-4:seq_len-1])
                lst_start.append(sequence[0:3])
                if start_position.span() != None:
                    if start_position.span() == (0, 3):
                        start_correct+=1
                    else:
                        print("start,{},{}".format(start_position.span(),h))
                        start_wrong += 1
                else:
                    continue
            #print(lst_start)
            print(set(lst_start))
            #print(lst_stop)
            print(set(lst_stop))
    #print(start_wrong)
        print(total)
        print(start_correct)
if __name__ == '__main__':
    main()
