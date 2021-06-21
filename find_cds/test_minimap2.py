from cutter_mpi import minimap2, convert_fasta
from glob import glob

ref_gene_files = glob('corrected_ref_cds/NC_006273/NC_006273_*')

# load query genomes
queries = convert_fasta('/home/sareh/data/pruned_genome/Pruned_nuc_NC_006273')


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
    rgene = convert_fasta(rgf)[0][1]
    for qh, qs in queries:
        qgene = minimap2(query=qs, refseq=rgene)
        p = pdist(qgene, rgene)
        print('{},{:1.3f}'.format(qh, 100*p))
    break
