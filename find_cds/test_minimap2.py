from cutter_mpi import minimap2, convert_fasta
from glob import glob
from cutter_mpi import mafft

ref_gene_files = glob('corrected_ref_cds/NC_000858/NC_00858_*')
#ref_gene_files = glob('corrected_ref_cds/NC_001414/NC_001414_*')
# runs with : corrected_ref_cds/NC_006273/NC_006273_*'
#all scores 0 : 'corrected_ref_cds/NC_001348/NC_001348_*'

# load query genomes
handle= open('/home/sareh/data/pruned_genome/Pruned_nuc_NC_000858')
#handle= open('/home/sareh/data/pruned_genome/Pruned_nuc_NC_001414')
# runs with : '/home/sareh/data/pruned_genome/Pruned_nuc_NC_006273'
#all scores 0: '/home/sareh/data/pruned_genome/Pruned_nuc_NC_001348'
queries = convert_fasta(handle)


def pdist(s1, s2):
    """
    Calculate the p-distance (proportion of nucleotide differences between
    two sequences)
    """
    ndiff = 0
    for i, nt1 in enumerate(s1):
        if nt1 in 'N-':
            continue
        nt2 = s2[i]
        if nt1 != nt2 and nt2 not in 'N-':
            ndiff += 1
    return ndiff/len(s1)


for rgf in ref_gene_files:
    handle = open(rgf)
    rgene = convert_fasta(handle)[0][1]
    for qh, qs in queries:
        #print(qh)
        qgene = minimap2(query=qs, refseq=rgene)
        qgene = mafft(query=qgene, ref=rgene, trim=False)
        ndiff = 0
        #print(rgene)
        #print(len(rgene))
        #print(qgene)
        #print(len(qgene))
        for i, nt1 in enumerate(qgene):
            #print("i is {}".format(i))
            #print("nt1 is {}".format(nt1))
            if nt1 in "N-":
                continue
            #if i < len(rgene):
            try:
                nt2 = rgene[i]
                #print("rgene[i] is {}".format(nt2))
                if nt1 != nt2 and nt2 not in "N-":
                    ndiff += 1
                #print(len(qgene)-len(rgene))
            except Exception as e:
                #print("i {}, nt1 {},nt2 {},qh {} \n {} \n {}".format(i,nt1,nt2,qh,qgene,rgene))
                print("i {}, nt1 {},nt2 {},qh {}".format(i,nt1,nt2,qh))
                print("qgene length:{}".format(len(qgene)))
                print("rgene length:{}".format(len(rgene)))
                print("ERROR : "+str(e))
        p = ndiff/len(qgene)
        #print(p)
        #p = pdist(qgene, rgene)
        print('{},{:1.3f}'.format(qh, 100*p))
        #NC_001348,KF853233.1,Human alphaherpesvirus 3,0.000
        #break
    #break
