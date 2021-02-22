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




"""
# read each file in Data/Pruned_FastTree (NC_0123456.tre)
# output a file for each file read with list of accn
import os
import re

def tre_accn_parser(file):
    #parses a .tre file to retrieve all accn numbers

    accns = []
    with open(file, 'r') as f:
        for line in f:
            x = re.findall('[A-Z]{2}[1-9]{6}.[1-9]', line)
            y = re.findall('[A-Z]{1}[1-9]{5}.[1-9]', line)
            # accn 1letter+5digits or 2letters+6digits
            # the y is only retriving one sequence (U55362.2)
        accns += x
        accns += y
    return (accns)

    # either read the entire file all at one and return findall
    # or make an empty list to append the content of findall for each line


def main():
    directory = '/Users/sarehchimeh/Data/Pruned_FastTree'
    for filename in os.listdir(directory):
        if filename.startswith('NC_'):
            # to avoid the ".DS_Store" files
            accn = filename[:9]
            print(accn)
            file_path = '/Users/sarehchimeh/Data/Pruned_FastTree/{}'.format(filename)
            with open('/Users/sarehchimeh/Data/pruned_accn/pruned_accn_{}'.format(accn), 'w') as outfile:
                values = tre_accn_parser(file_path)
                outfile.write(values)


if __name__ == "__main__":
    main()
"""



"""
#check all lists should be less than 100 elements
dir= '/Users/sarehchimeh/Data/pruned_accn'
for filename in os.listdir(dir):
    if filename.startswith('NC_'):
        path = '{}/{}'.format(dir,filename)
        with open(path, 'r') as f:
            for line in f:
                if len(line) >100:
                    print(filename)
 """

#write a bash script iterate through the files and wc it with an if file to check their<100