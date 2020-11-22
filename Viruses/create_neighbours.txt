import csv

Neighbours = '/Users/sarehchimeh/Data/NCBI_Neighbours.txt'

def main():
    with open(Neighbours,'r') as file:
        title = file.readline()
        line_two = file.readline()
        header = line_two[12:]
        with open('/Users/sarehchimeh/Data/Neighbours.txt','w') as of:
            of.write(title)
            of.write(header)
            for line in file:
                of.write(line)

if __name__ == "__main__":
    main()
