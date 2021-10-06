from mpi4py import MPI
import sys
import math
import os
from itertools import permutations

path = "/home/sareh/surfaces/fubar/json_names.txt"

# MPI
my_number = MPI.COMM_WORLD.Get_rank() # ID of this process
total_number = MPI.COMM_WORLD.Get_size() # total number of processes

def file_names(path):
    """
    takes a file path and return each line as a list item
    """
    name_lst = []
    file = open(path,'r')
    for line in file:
        name_lst.append(line.strip("\n"))
    file.close()
    return(name_lst)

def main():
    count = 0
    # get a list of all names of gene files
    name_lst = file_names(path)
    # get all permuations of the pair of genes (order matters)
    perm = permutations(name_lst, 2)

    for pair in list(perm):
        count += 1 
        if count % total_number != my_number:
            continue
        run = ("Rscript json_pearsoncor.r {} {}".format(pair[0],pair[1]))
        os.system(run)
        # Rscript json_pearsoncor.r NC_000858_NP_049558.1.clean.FUBAR.json NC_000858_NP_049559.1.clean.FUBAR.json

if __name__ == '__main__':
    main()


"""
    # mpi run
    for i in range(1, 101):
        if i % nprocs == my_rank:
            os.system(run)
"""
