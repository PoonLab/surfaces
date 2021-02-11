Last login: Thu Feb 11 10:27:07 on ttys003
(base) sarehchimeh@Sarehs-MacBook-Air ~ % cd 
(base) sarehchimeh@Sarehs-MacBook-Air ~ % cd PycharmProjects 
(base) sarehchimeh@Sarehs-MacBook-Air PycharmProjects % ls
Surfaces	git		pythonProject	pythonProject1
(base) sarehchimeh@Sarehs-MacBook-Air PycharmProjects % cd Surfaces 
(base) sarehchimeh@Sarehs-MacBook-Air Surfaces % ls
Downsizing		Viruses			practice.py
README.md		overlapping_genes	retrieve_sequences
(base) sarehchimeh@Sarehs-MacBook-Air Surfaces % git remote -v 
origin	https://github.com/Sarehchimeh/Surfaces.git (fetch)
origin	https://github.com/Sarehchimeh/Surfaces.git (push)
(base) sarehchimeh@Sarehs-MacBook-Air Surfaces % git remote set-url poonlab https://github.com/PoonLab/Surfaces.git
fatal: No such remote 'poonlab'
(base) sarehchimeh@Sarehs-MacBook-Air Surfaces % git remote add poonlab https://github.com/PoonLab/Surfaces.git
(base) sarehchimeh@Sarehs-MacBook-Air Surfaces % git remote -v
origin	https://github.com/Sarehchimeh/Surfaces.git (fetch)
origin	https://github.com/Sarehchimeh/Surfaces.git (push)
poonlab	https://github.com/PoonLab/Surfaces.git (fetch)
poonlab	https://github.com/PoonLab/Surfaces.git (push)
(base) sarehchimeh@Sarehs-MacBook-Air Surfaces % git push poonlab master
error: src refspec master does not match any
error: failed to push some refs to 'https://github.com/PoonLab/Surfaces.git'
(base) sarehchimeh@Sarehs-MacBook-Air Surfaces % ssh sareh@129.100.26.213 
sareh@129.100.26.213's password: 
Permission denied, please try again.
sareh@129.100.26.213's password: 
Welcome to Ubuntu 18.04.5 LTS (GNU/Linux 4.15.0-132-generic x86_64)

 * Documentation:  https://help.ubuntu.com
 * Management:     https://landscape.canonical.com
 * Support:        https://ubuntu.com/advantage

 * Canonical Livepatch is available for installation.
   - Reduce system reboots and improve kernel security. Activate at:
     https://ubuntu.com/livepatch

58 packages can be updated.
3 updates are security updates.

New release '20.04.2 LTS' available.
Run 'do-release-upgrade' to upgrade to it.

                              

#! /usr/bin/python3
#prune a fasta file from (Data/CDS nucleotide sequence)
# to keep only the sequence entries of the accn sequences in Data/pruned_accn
# input Data/CDS nucleotide sequence & Data/pruned_accn directories

import os
import re
from Bio import SeqIO

def get_accn(file):
    """
    Stream through txt file retrieve accns into a list
    """
    lst_accns = []
    with open(file, 'r') as f:
        for line in f:
            line = line.strip('\n')
            lst_accns.append(line)
    a = tuple(lst_accns)
    return(a)

# if "accn" in tuple:

def main():
    accn_directory = '/home/sareh/data/pruned_accns'
    cds_directory = '/home/sareh/data/CDS_nucleotide_seq'
    for filename1 in os.listdir(accn_directory):
        print('filename1 is '+ filename1)
        #filename1='pruned_accn_NC_007605'
        accn1=filename1[12:]
        print('accn1 ' + accn1)
        path1 = '{}/{}'.format(accn_directory, filename1)
        path2 = '{}/{}.txt'.format(cds_directory, accn1)
        print(path2)
        # filename2='NC_001526.fasta'
        accessions= get_accn(path1)
        print(accessions)
        #all accessions in a tuple
        with open('/home/sareh/data/Pruned_CDS/Pruned_CDS_{}'.format(accn1),'w') as o$
            pruned_record_count = 0
            record_count = 0
            with open(path2,'r') as f2:
                for record in SeqIO.parse(path2, 'fasta'):
                    record_count += 1
                    id = []
                    a = record.id
                    x = re.findall('[A-Z]{1,3}[0-9]{5,6}.[0-9]', a)
                    # accns in a list
                    for i in x:
                        if i in accessions:
                            pruned_record_count +=1
                            SeqIO.write(record,out_file,'fasta')
        print('record count: '+ str(record_count))
        print('pruned record count: '+ str(pruned_record_count))
                                    # write the next line into a file (record sequenc$

if __name__ == "__main__":
    main()



