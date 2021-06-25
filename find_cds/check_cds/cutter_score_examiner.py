"""
cutter_scores:
header,align.score
NC_000858,LC210029.1,Human T-cell leukemia virus type I\,4.565891472868217
writes out a csv file with a row summarizing each file (gene):
("virus,max,min,zero,one,two,three,four\n")
"""

import os
import re
import csv

directory= "/home/sareh/surfaces/find_cds/cutter_scores"
outfile_path = '/home/sareh/surfaces/find_cds/check_cds/examine_cutter_scores.csv'

#directory= "/home/sareh/surfaces/find_cds/check_cds/test/backup_cutter_scores"
#outfile_path = '/home/sareh/surfaces/find_cds/check_cds/test/examine_cutter_scores.csv'

#print("virus,max,min,zero,one,two,three,four")

num_files = 0
num_alignments = 0
all_zero = 0
all_one = 0
all_two = 0

outfile = open(outfile_path,'w')
#Header of csv file
outfile.write("virus,max,min,zero,one,two,three,four\n")

for file in os.listdir(directory):
    #print(file)
    num_files += 1
    #reading each scores file for each gene
    path =("{}/{}".format(directory,file))
    with open(path, "r") as f:
        total = 0
        row_count = 0
        four = 0
        three = 0
        two = 0
        one = 0
        zero = 0
        list = []
        row_count = 0
        csv_reader = csv.reader(f,delimiter=',')
        next(csv_reader,None)
        for row in csv_reader:
            #print(row)
            num_alignments += 1
            s = float(row[3])
            row_count += 1
            list.append(s)
            if s >= 4:
                four += 1
            elif s >= 3:
                three += 1
            elif s >= 2:
                two += 1
                all_two +=1
            elif s >= 1 and s<2:
                one += 1
                all_one +=1
            elif s < 1:
                print(row)
                zero += 1
                all_zero += 1
         
       
            #print("zero is {}".format(zero))
            #print("one is {}".format(one))
            #print("two is {}".format(two))
            #print("three is {}".format(three))
            #print("four is {}".format(four))
        if len(list) > 0:
            outfile.write("{},{:.1f},{:.2f},{},{},{},{},{}\n".format(file,max(list),min(list),zero,one,two,three,four))
        #NC_007605_YP_401637.1,5.0,0.59,31,0,0,0,69
            
print("num of files is {}".format(num_files))
print("num of alignments is {}".format(num_alignments))
print("all zero {}".format(all_zero))
print("all one {}".format(all_one))
print("all two {}".format(all_two))

outfile.close
