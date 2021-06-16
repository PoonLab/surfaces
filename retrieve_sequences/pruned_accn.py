#! /usr/bin/python
#read pruned.tre to get the list of accns after pruning
#read each file in Data/Pruned_FastTree (NC_0123456.tre)
#output a file for each file read with list of accn
import os
import re

directory= '/home/sareh/data/Pruned_FastTree'
for filename in os.listdir(directory):
    if filename.startswith('NC_'):
           # to avoid the ".DS_Store" files
           #filename "NC_004718_pruned.tre"
        accn=filename[:9]
        file_path = '{}/{}'.format(directory,filename)
        with open('/home/sareh/data/pruned_accns/pruned_accn_{}'.format(accn), 'w') as outfile:
            accns = []
            print(accn)
            with open(file_path, 'r') as f:
                for line in f:
                    print(line)
                    #why isn't this printing
                    x = re.findall('[A-Z]{1,3}[0-9]{5,6}.[0-9]', line)
                        # accn 2letters+6digits.digit
                        # accn 1letter+5digits.digit
                        # accn 3letters 
                        # the y is only retriving one sequence (U55362.2)
                accns += x
                #accns is a list here
            for item in accns:
                outfile.write('{}\n'.format(item))
            #write each accn in a line

