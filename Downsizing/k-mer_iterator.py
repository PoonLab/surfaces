#k-mer iterator
#python3 k-mer.py /home/sareh/data/Pruned_CDS/Pruned_CDS_NC_005148 Kmer_outfle_NC_005148.csv
import os
count = 0
file='/home/sareh/data/Pruned_CDS'
with open ('/home/sareh/data/k-mer_iterator','w') as f:
    for filename in os.listdir(file):
        if 'NC_' in filename:
            print(filename)
            accn=filename[11:]
            print(accn)
            outfile= accn +'_K-mer.csv'
            print(outfile)
            f.write('python3 k-mer.py /home/sareh/data/Pruned_CDS/Pruned_CDS_{} /home/sareh/data/K-mer_output/{}\n'.format(accn, outfile))
            count+=1
print(count)
