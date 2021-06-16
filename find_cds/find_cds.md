# find_cds

## cutter.py 
Cutting out each cds/gene from neighbour seqeucnes (pruned, 100/virus species) from the ref sequence's cds annotation in NCBI
cutter_cds

## Why
The NCBI cds annotations was not consistant enough to retrieve all cds for the pruned seqeunces 

## Goal
To get a fasta file for each of the genes in each virus species  
(100 entries in each fasta file)
Use codon aware alignment to enter gaps 
1. Convert to aa | nuc_to_aa.py | translate_nuc function from gotoh2 Utils 
2. align aa | mafft
aa_cds 
3. insert codon aware gaps | apply_prot_tonuc from gotoh2 Utils
gap_nuc_cds

## Issues 
* Some cds used in alignemnt for the cutter script are antisense (-1) 
	if the header of nuc_refseq_cds is (-1) then take the reverse complement and run cutter 
	1. find_rev.py -> rev_nuc_refseq_cds 
	   198 files are (-1)
 
* Cutter script takes a week to run so need cluster computing   
