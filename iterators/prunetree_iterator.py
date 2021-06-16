#! /usr/bin/bash
#Prunetree iterator
#python3 prunetree.py NC_000858.tre 100 > prune_test.txt

import os
count = 0
file='/home/sareh/data/FastTree_output'
with open ('/home/sareh/data/prunetree_iterator.sh','w') as f:
    for filename in os.listdir(file):
        if 'NC_' in filename:
   #filename=FastTree_outputNC_044047.tre
            accn=filename[15:-4]
            print(accn)
            outfile= accn +'_pruned.tre'
            f.write('python3 /home/sareh/script/prunetree.py /home/sareh/data/FastTree_output/{} 100 > /home/sareh/data/Pruned_FastTree/{}\n'.format(filename, outfile))
            count+=1
print(count)

