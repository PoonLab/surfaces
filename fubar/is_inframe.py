dir = "/home/sareh/re-do/100/100cds/cut_cds/"
#dir = "/home/sareh/hiv/FUBAR/hyphy_clean" # clean.NC_001802_NP_057857.2

import Bio
from Bio import AlignIO
from Bio.Align import AlignInfo
from Bio import Entrez, SeqIO, SeqFeature
from csv import DictWriter
from time import sleep
import sys
import os
import re
# Email to be identified in NCBI
Entrez.email = 'sbagher4@uwo.ca'

codon_dict = {'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L',
                'TCT':'S', 'TCC':'S', 'TCA':'S', 'TCG':'S',
                'TAT':'Y', 'TAC':'Y', 'TAA':'*', 'TAG':'*',
                'TGT':'C', 'TGC':'C', 'TGA':'*', 'TGG':'W',
                'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L',
                'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P',
                'CAT':'H', 'CAC':'H', 'CAA':'Q', 'CAG':'Q', 
                'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R',
                'ATT':'I', 'ATC':'I', 'ATA':'I', 'ATG':'M',
                'ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T', 
                'AAT':'N', 'AAC':'N', 'AAA':'K', 'AAG':'K',
                'AGT':'S', 'AGC':'S', 'AGA':'R', 'AGG':'R',
                'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V',
                'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',
                'GAT':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E',
                'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G',
                '---':'-', 'XXX':'?'}

mixture_regex = re.compile('[WRKYSMBDHVN-]')

mixture_dict = {'W':'AT', 'R':'AG', 'K':'GT', 'Y':'CT', 'S':'CG', 
                'M':'AC', 'V':'AGC', 'H':'ATC', 'D':'ATG', 
                'B':'TGC', 'N':'ATGC', '-':'ATGC'}

#mixture_dict_2 =  [ (set(v), k) for k, v in mixture_dict.iteritems() ]
ambig_dict = dict(("".join(sorted(v)), k) for k, v in mixture_dict.items())

def translate_nuc (seq, offset, resolve=False, return_list=False):
	"""
	Translate nucleotide sequence into amino acid sequence.
		offset by X shifts sequence to the right by X bases
	Synonymous nucleotide mixtures are resolved to the corresponding residue.
	Nonsynonymous nucleotide mixtures are encoded with '?' 
	"""
	
	seq = '-'*offset + seq
	
	aa_list = []
	aa_seq = ''	# use to align against reference, for resolving indels
	
	# loop over codon sites in nucleotide sequence
	for codon_site in range(0, len(seq), 3):
                
		codon = seq[codon_site:codon_site+3]
		
		if len(codon) < 3:
			break
		
		# note that we're willing to handle a single missing nucleotide as an ambiguity
		if codon.count('-') > 1 or '?' in codon:
			if codon == '---':	# don't bother to translate incomplete codons
				aa_seq += '-'
				aa_list.append(['-'])
			else:
				aa_seq += '?'
				aa_list.append(['?'])
			continue
		
		# look for nucleotide mixtures in codon, resolve to alternative codons if found
		num_mixtures = len(mixture_regex.findall(codon))
		
		if num_mixtures == 0:
			aa = codon_dict[codon]
			aa_seq += aa
			aa_list.append([aa])
			
		elif num_mixtures == 1:
			resolved_AAs = []
			for pos in range(3):
				if codon[pos] in mixture_dict.keys():
					for r in mixture_dict[codon[pos]]:
						rcodon = codon[0:pos] + r + codon[(pos+1):]
						if codon_dict[rcodon] not in resolved_AAs:
							resolved_AAs.append(codon_dict[rcodon])
							
			aa_list.append(resolved_AAs)
			
			if len(resolved_AAs) > 1:
				if resolve:
					# for purposes of aligning AA sequences
					# it is better to have one of the resolutions
					# than a completely ambiguous '?'
					aa_seq += resolved_AAs[0]
				else:
					aa_seq += '?'
			else:
				aa_seq += resolved_AAs[0]
				
		else:
			aa_seq += '?'
			aa_list.append(['?'])
			
	if return_list:
		return aa_list

	return aa_seq

def translate(seq):
    table = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
        'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W',
    }
    protein = ""
    if len(seq) % 3 == 0:
        for i in range(0, len(seq), 3):
            try:
                codon = seq[i:i + 3]
                protein += table[codon]
            except KeyError:
                protein += "X"
    return protein


def pdist(s1, s2):
    """
    Calculate the p-distance (proportion of nucleotide differences between
    two sequences)
    """
    ndiff = 0
    diff_site =[]
    for i, nt1 in enumerate(s1):
        if nt1 in 'NX':
            continue
        try:
            nt2 = s2[i]
            if nt1 != nt2:
                ndiff += 1
                diff_site.append(i)
        except:
            continue
    #return ndiff / len(s1)
    return diff_site

def retrieve_gid(accn):
    """
    Find records according to genbank accession number
    """
    handle = Entrez.esearch(db='protein', term=accn, retmax=1)
    response = Entrez.read(handle)
    if response['Count'] == '0':
        # retry query
        handle = Entrez.esearch(db='protein', term=accn, retmax=1)
        response = Entrez.read(handle)
        if response['Count'] == '0':
            return None

    return response['IdList'][0]

def retrieve_record(gid):
    handle = Entrez.efetch(db='protein', rettype='gb', retmode='text', id=gid)
    gb = SeqIO.parse(handle, format='genbank')
    return next(gb)


def main():
    outfile_cons_path = ("{}consensus0.5.fa".format(dir))
    outfile_csv_path = ("{}is_inframe0.5.csv".format(dir))
    outfile_nuc_path = ("{}consensus0.5.nuc.fa".format(dir))
    #print(outfile_cons_path)
    #print(outfile_csv_path)

    outfile_cons = open(outfile_cons_path, "w")
    outfile_csv = open(outfile_csv_path, "w")
    outfile_nuc = open(outfile_nuc_path,"w")    

    outfile_csv.write( "{},{},{},{},{},{},{}\n".format("id", "len_consensus", "len_consensusaa",
                                         "len_aa","num_X","percentX","hamming_dist",
                                         "list_diff", "consecutive"))

    for file in os.listdir(dir):
        path = os.path.join(dir, file)
        print(file)

        #id from file name
        b = file.replace("clean.","") #begining of the file name
        b = b.replace(".fasta", "") #end of the file name
        c = re.match("[A-Z]{2}\_[0-9]{6}\_", b) #the virus accn
        virus = c.group()
        ID = b.replace(virus,"")
        #print(ID)

        # Opening the fasta msa file
        align = AlignIO.read(path, "fasta")
        # Instantiate SummaryInfo class and pass 'align' to it.
        summary_align = AlignInfo.SummaryInfo(align)
        consensus = summary_align.dumb_consensus(threshold=0.5)
        #print(consensus)
        #print(len(consensus))
        cons_aa = translate(consensus)        
        #print(cons_aa)

        # Retrieve protein sequence
        gid = retrieve_gid(ID)
        if gid is None:
            print('Warning, failed to retrieve gid for {}'.format(ID))
        else:
            record = retrieve_record(gid)
            aaseq = record.seq

        # P-dist calculation
        pdistance= pdist(cons_aa,aaseq)
        consecutive = "FALSE"
        if len(pdistance) == 0:
            hamming_dist = 0
            p_dist = 0
        else:
            hamming_dist = len(pdistance)
            p_dist = hamming_dist/len(cons_aa)
            for i in range(0,len(pdistance)-3):
                if (pdistance[i+1])-(pdistance[i])==1 & (pdistance[i+2])-(pdistance[i+1])==1 & (pdistance[i+3])-(pdistance[i+2])==1:
                    #print("its consecutive")
                    consecutive = "TRUE"
        perctX = cons_aa.count("X")/len(aaseq)

        #print(file + "," + str(cons_aa.count("X")) + "," + str(len(aaseq)) + "," + str(perct))
        #print(">{}\n{}\n>{}\n{}\n".format(ID, cons_aa,ID,aaseq))

        # Write out to a file
        outfile_nuc.write(">{}\n{}\n".format(ID, consensus))
        outfile_cons.write(">{}\n{}\n>{}\n{}\n".format(ID, cons_aa,ID,aaseq))
        outfile_csv.write("{},{},{},{},{},{},{},{}\n".format(ID,len(consensus),len(cons_aa),
                                                     len(aaseq),cons_aa.count("X"),perctX,hamming_dist,
                                                     p_dist,consecutive))

    outfile_cons.close()
    outfile_csv.close()
    outfile_nuc.close()


if __name__ == "__main__":
    main()


