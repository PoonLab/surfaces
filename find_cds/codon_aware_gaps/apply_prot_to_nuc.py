"""
takes gaps in the aligned aa sequences and applies that to the nucloetide sequences 
aa_directory of aa aligned fasta files 
nuc_directory of unaligned/original nuc fasta files 
"""

import re
import os

aa_directory = "/home/sareh/surfaces/find_cds/codon_aware_gaps/aligned_cds"
nuc_directory = "/home/sareh/surfaces/find_cds/codon_aware_gaps/cut_cds"
out_directory = "/home/sareh/surfaces/find_cds/codon_aware_gaps/gap_nuc_CDS"
    
def parse_fasta (handle):
    """
    Parse open file as FASTA, return dictionary of 
    headers and sequences as key-value pairs.
    """
    res = {}
    sequence = ''
    for i in handle:
        if i[0] == '$': # skip h info
            continue
        elif i[0] == '>' or i[0] == '#':
            if len(sequence) > 0:
                res.update({h: sequence})
                sequence = ''   # reset containers
                h = i.strip('\n')[1:]
            else:
                h = i.strip('\n')[1:]
        else:
            sequence += i.strip('\n').upper()
    res.update({h: sequence})
    return res


def apply_prot_to_nuc(aligned_prot, nuc):
    """
    Use aligned protein sequence to update the corresponding nucleotide
    sequence.
    :param aligned_prot:  str, aligned amino acid sequence
    :param nuc:  str, original nucleotide sequence
    :return: str, nucleotide sequence with codon gaps
    """
    res = ''
    i = 0
    for aa in aligned_prot:
        if aa == '-':
            res += '---'
            continue
        res += nuc[i:(i + 3)]
        i += 3
    return res


def main():
    file_lst=[]
    for nuc_file in os.listdir(nuc_directory):
        print(nuc_file) #cutter_NC_003045_NP_150080.1
        name = re.sub("cutter_","",nuc_file)
        file_lst.append(name) #print(name)
    #print(file_lst) #['NC_003045_NP_150080.1', 'NC_003045_NP_150079.1']

    for file in file_lst:

        #print(file) #NC_003045_NP_150080.1
        aa_path_name = "aaCDS_cutter_{}".format(file)
        aa_path = os.path.join(aa_directory,aa_path_name)
        nuc_path_name = "cutter_{}".format(file)
        nuc_path = os.path.join(nuc_directory,nuc_path_name)

        try:
            with open (aa_path,"r") as aa_f:
                aa_dict = parse_fasta(aa_f)
                #print(len(aa_dict))
        except IOError:
            print("file {} does not exsist".format(aa_path))

        try:
            with open (nuc_path,'r') as nuc_f:
                nuc_dict = parse_fasta(nuc_f)
                #print(len(nuc_dict))
        except IOError:
            print("file {} does not exsist".format(nuc_path))
        
        out_path = os.path.join(out_directory,file)
        with open (out_path,'w+') as outfile:
            for key in aa_dict:
                aligned_aa = aa_dict[key]
                nuc = nuc_dict[key]
                new_seq=apply_prot_to_nuc(aligned_aa,nuc)
                outfile.write("{}\n{}\n".format(key,new_seq))

if __name__ == '__main__':
    main()
