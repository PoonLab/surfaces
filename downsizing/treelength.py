"""
Read in Newick files and calculate the tree length 
treelength: sum of branch lengths in the tree
"""

import os 
import re

#dir = "/home/sareh/data/fasttree_output"
dir = "/home/sareh/surfaces/find_cds/working/tree_for_fasttree"
outfile_path = "/home/sareh/data/neighbor_cds/treelength_cds100.csv"
outfile = open(outfile_path,"w")

pattern = re.compile("^[0-9|.]*")

for file_name in os.listdir(dir):
    #virus_name=file_name.replace("FastTree_output","")
    #virus_name=virus_name.replace(".tre","")
    name = file_name.replace("cutter_cds_","")
    #print(file_name) #FastTree_outputNC_038889.tre
    #print(virus_name) #NC_038889
    #file_path ="/home/sareh/data/fasttree_output/FastTree_outputNC_038889.tre"
    file_path =os.path.join(dir,file_name)
    file = open(file_path,"r")
    t = file.readline()

    lengths=[]
    for l in t.split(":"):
        branch_len = pattern.findall(l)
        int_list = [lengths.append(i) for i in branch_len]

    lengths = list(filter(None, lengths)) # Taking the empty items out of list 
    lengths = list(map(float, lengths)) # Turning str to numbers

    outfile.write("{},{}\n".format(name,sum(lengths)))
    
