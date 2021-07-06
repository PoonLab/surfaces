#combine testminimap and cutter_mpi.py
from glob import glob

#mpi
from mpi4py import MPI
the_world = MPI.COMM_WORLD
my_number = the_world.Get_rank()
total_number = the_world.Get_size()

def get_accn(file):
    """
    takes in file path
    Stream through txt file retrieve accns into a list
    """
    lst_accns = []
    with open(file, 'r') as f:
        for line in f:
            line = line.strip('\n')
            lst_accns.append(line)
    a = tuple(lst_accns)
    return(a)

def convert_fasta (file):
    """
    takes in file path 
    """
    result = []
    sequence = ''
    with open(file,'r') as handle:
        for line in handle:
            if line.startswith('$'): # skip header line
                continue
            elif line.startswith('>') or line.startswith('#'):
                if len(sequence) > 0:
                    result.append([h,sequence])
                    sequence = ''   # reset
                h = line.strip('>#\n')
            else:
                sequence += line.strip('\n').upper()
            
        result.append([h,sequence]) # handle last entry
        return result

def pdist(s1, s2):
    """
    Calculate the p-distance (proportion of nucleotide differences between
    two sequences)
    """
    ndiff = 0
    for i, nt1 in enumerate(s1):
        if nt1 in 'N-':
            continue
        nt2 = s2[i]
        if nt1 != nt2 and nt2 not in 'N-':
            ndiff += 1
    return ndiff/len(s1)
    
#Directories 
ref_home_directory = '/home/sareh/surfaces/find_cds/corrected_ref_cds'
query_directory = '/home/sareh/data/pruned_genome'
accession_file = '/home/sareh/virus_corrected_ref_cds.txt'

def main():
    accns = get_accn(accession_file) #all of the files to examine
    for accn in accns:
        ref_gene_files = glob('/home/sareh/surfaces/find_cds/corrected_ref_cds/{}/{}_*'.format(accn,accn))
        query_file_path = ('/home/sareh/data/pruned_genome/Pruned_nuc_{}'.format(accn))
        #query_file_handle = open(query_file_path)
        queries = convert_fasta(query_file_handle)
        
        #Each gene files
        for rgf in ref_gene_files:
            #rfg_handle = open(rfg)
            rgene = convert_fasta(rgf)[0][1] #the reference cds seqeunce
            
            #Each query header and sequence
            for qh, qs in queries:
                #print(qh)
                qgene = minimap2(query=qs, refseq=rgene) #each alignment 
                ndiff = 0
                #print(rgene) #print(len(rgene)) 
                #print(qgene) #print(len(qgene))
                
                try:
                    p = pdist(qgene, rgene)
                    
                except Exception as e:
                    print("i {}, nt1 {},nt2 {},qh {} \n {} \n {}".format(i,nt1,nt2,qh,qgene,rgene))
                    print("qgene length:{}".format(len(qgene)))
                    print("rgene length:{}".format(len(rgene)))
                    print("ERROR : "+str(e))
                
                """
                #Alignment Score Calculation (P-distance)
                for i, nt1 in enumerate(qgene):
                    #print("i is {}".format(i))  #print("nt1 is {}".format(nt1))
                    if nt1 in "N-":
                        continue
                    #if i < len(rgene):
                    try:
                        nt2 = rgene[i]
                        #print("rgene[i] is {}".format(nt2))
                        if nt1 != nt2 and nt2 not in "N-":
                            ndiff += 1
                 
                    #Print Errors 
                    except Exception as e:
                        print("i {}, nt1 {},nt2 {},qh {} \n {} \n {}".format(i,nt1,nt2,qh,qgene,rgene))
                        print("qgene length:{}".format(len(qgene)))
                        print("rgene length:{}".format(len(rgene)))
                        print("ERROR : "+str(e))
                        #print(len(qgene)-len(rgene))
                p = ndiff/len(qgene)
                """ 
                
                print('{},{:1.3f}'.format(qh, 100*p))
                
            #break
            
        break
  
if __name__ == "__main__":
    main()
