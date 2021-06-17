#cutter_scores
import os
import re

directory= "/home/sareh/cutter/cutter_scores"

#print("virus,max,min,zero,one,two,three,four")

num_files = 0
num_alignments = 0
empty_files = 0 
all_zero = 0
all_one = 0
all_two = 0
outfile = open('/home/sareh/cutter/examine_cutter_new.csv','w')
outfile.write("virus,max,min,zero,one,two,three,four\n")
for file in os.listdir(directory):
    #print(file)
    num_files += 1
    path =("{}/{}".format(directory,file))
    with open(path, "r") as f:
        total = 0
        four = 0
        three = 0
        two = 0
        one = 0
        zero = 0
        list = []
        for line in f:
            if re.match("^[a-zA-Z]", line):
                continue 
            else:
                total += 1
                line = (line.strip("\n"))
                list.append(float(line))
        #print(list)

        if total ==0:
            print("{} is empty".format(file))
            empty_files += 1
            continue
        else:
            #print(file)
            #print("max is {}".format(max(list)))
            #print("min is {}".format(min(list)))
            for i in list:
                num_alignments +=1
                if i >= 4:
                    four += 1
                elif i >= 3:
                    three += 1
                elif i >= 2:
                    two += 1
                    all_two +=1
                elif i >= 1 and i<2:
                    one += 1
                    all_one +=1
                elif i < 1:
                    print(i)
                    zero += 1
                    all_zero += 1

            #print("zero is {}".format(zero))
            #print("one is {}".format(one))
            #print("two is {}".format(two))
            #print("three is {}".format(three))
            #print("four is {}".format(four))

            outfile.write("{},{:.1f},{:.2f},{},{},{},{},{}\n".format(file,max(list),min(list),zero,one,two,three,four))
            #NC_007605_YP_401637.1,5.0,0.59,31,0,0,0,69

print("num of files is {}".format(num_files))
print("num of alignments is {}".format(num_alignments))
print("all zero {}".format(all_zero))
print("all one {}".format(all_one))
print("all two {}".format(all_two))
print("empty files {}".format(empty_files))
