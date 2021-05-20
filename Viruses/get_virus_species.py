"""
Parse the NCBI virus table from : https://www.ncbi.nlm.nih.gov/genomes/GenomesGroup.cgi?taxid=10239&host=human
Infile: download 1. "the table" 2. "Genome Neighbors Data" (filter for human viruses)
outfiles:
ref_accns.txt - retrieve the Ref (Representative) accessions with >= 3 proteins and >= 100 neighbours available
neighbor_accns.txt - retrieve all the neighbour accessions for the viruses with  >= 3 proteins and >= 100 neighbours available
*the header in the neighbour file don't show up properly sometimes! (delete ##column)
"""
#1. "the table"
path_table = "/Users/sareh/Desktop/ncbi_table.tbl"
#2. "Genome Neighbors Data"
path_neighbors = "/Users/sareh/Desktop/ncbi_neighbors.nbr"

ref_outfile_path = "/Users/sareh/Desktop/ref_accns.txt"
neighbors_outfile_path = "/Users/sareh/Desktop/neighbor_accns.txt"
neighbors_directory_path = "/Users/sareh/Desktop/Neighbors"

import csv
import os

ref_list = []

def parse_table(table):
    """
    :param table: file from ncbi 1. "the table"
    :return: list of accessions that meet the criteria
    """
    with open (path_table, 'r') as table:
        title = table.readline()
        colnames = table.readline()
        table_csv = csv.reader(table, delimiter='\t')
        count_row = 0
        count_virus = 0
        for row in table_csv:
            count_row +=  1
            #remove entries begining in tabs (plus sign on the website)
            if len(row)<2 or row[7]=="-" or row[6]=="-":
                continue
            #remove the rows with missing neighbour number or Representative accens
            elif int(row[7]) >= 100 and int(row[6]) >= 3:
                #print(row)
                count_virus +=1
                ref_list.append(row[1])

    return(ref_list,count_virus,count_row)


def main():
    count_ref = 0
    ref_list,count_virus,count_row=parse_table(path_table)
    print(ref_list)

    #Write out the ref accens
    with open(ref_outfile_path, 'w') as ref_outfile:
        for ref in ref_list:
            count_ref +=1
            ref_outfile.write(ref+"\n")

    count = 0
    multiple =[]
    with open(path_neighbors, 'r') as file:
        title = file.readline()
        #read into a dictionary the colnames are the keys to each dict
        reader = csv.DictReader(file, delimiter='\t')
        # opening the Neighbours file as a csv
        columns = ['Representative', 'Neighbor', 'Host', 'Selected lineage', 'Taxonomy name', 'Segment name']
        for virus in ref_list:
            filepath = '{}/{}_neighbors.txt'.format(neighbors_directory_path,virus)
            with open(neighbors_outfile_path, 'w+') as all_of:
                all_writer = csv.DictWriter(all_of, fieldnames=columns, delimiter='\t')
                with open(filepath, 'w') as ind_of:
                    #of.write(title)
                    ind_writer = csv.DictWriter(ind_of, fieldnames=columns, delimiter='\t')
                    ind_writer.writeheader()
                    file.seek(0)
                    next(file)
                    for row in reader:
                        #print(row)
                        #print(row)
                        print(row['Representative'])
                        print(row['Neighbor'])
                        rep_lst = row['Representative'].split(',')
                        print(rep_lst)
                        #print(rep_lst)
                        if virus in rep_lst:
                            ind_writer.writerow(row)
                            all_writer.writerow(row)
                            count += 1
    print("total rows {}".format(count_row))
    print("total_viruses {}".format(count_virus))
    print(count)
if __name__ == "__main__":
    main()
