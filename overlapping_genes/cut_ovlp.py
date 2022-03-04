"""
testing cutting ovlp regions
NC_001355,"NP_040297.1",529,1,"NP_040296.1",101,1,297,453,26,+2
"""
import os
import tempfile
import subprocess

def all_positions(loc):
    """
    loc =  2989:6643 OR 2989:6643;4566:6233
    """
    lst = loc.split(";")
    for i in lst:
        l = i.split(":")
        start = l[0]
        stop = l[1]
    return(start,stop)

def convert_fasta (handle):
    """
    takes in handel to an open fasta file
    """
    result = []
    sequence = ''
    # handle = open(file,'r')
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

def mafft(query, ref, trim=True):
    handle = tempfile.NamedTemporaryFile(delete=False)
    s = '>ref\n{}\n>query\n{}\n'.format(ref, query)
    handle.write(s.encode('utf-8'))
    handle.close()

    # call MAFFT on temporary file
    stdout = subprocess.check_output(['mafft', '--quiet', handle.name]) #Failed to align Error
    stdout = stdout.decode('utf-8')
    result = convert_fasta(stdout.split('\n'))
    aligned_ref = result[0][1]
    aligned_query = result[1][1]

    # trim aligned query sequence to extent of reference
    if trim:
        left, right = get_boundaries(aligned_ref)
        aligned_query = aligned_query[left:right]
    os.remove(handle.name)  # clean up
    return(aligned_query, aligned_ref)

def trimed_seq(org_seq, lst_ovlp,start):
    new_seq=""
    for i, nt  in enumerate(org_seq):
        genome_i = int(i) + int(start)
        if genome_i in lst_ovlp:
            continue
        else:
            new_seq += nt
    return(new_seq)

def main():
    removed = 0
    for virus in os.listdir(dir):
    #virus = "NC_001355"
        print(virus)
        virus_ovlp = []
        path = os.path.join(dir,virus)
        handle = open(path,"r")
        next(handle) #skip header line
        for line in handle:
            #"NC_001357,NP_040314.1,3417,1,NP_040313.1,2816,1,267,1098,268,+1"
            line_lst = line.split(",")

            gene1 = list(range(int(line_lst[2]) , int(line_lst[2]) + int(line_lst[7]) + 1 ))
            gene2 = list(range(int(line_lst[5]) , int(line_lst[5]) + int(line_lst[8]) + 1 ))
            line_ovlp = list(set(gene1).intersection(gene2))
            virus_ovlp.extend(line_ovlp)
        all_ovlp = list(set(virus_ovlp))
        print(len(all_ovlp))
        #exec(name + " = unq_lst")

        ref_dir = "/home/sareh/data/sequences/ref_cds"
        ref_virus_path = os.path.join(ref_dir, virus)
        if os.path.isdir(ref_virus_path):
            for ref_cds in os.listdir(ref_virus_path):
                print(ref_cds)
                ref_cds_path = os.path.join(ref_virus_path, ref_cds)
                ref_cds_handle = open(ref_cds_path, "r")

                header = ref_cds_handle.readline()
                h = ((header.strip(">")).strip("\n")).strip('"')
                info = h.split(",")

                start,stop = all_positions(info[4])
                cds_coord = list(range(int(start),int(stop)))

                print(len(cds_coord))
                cds_ovlp = set(all_ovlp).intersection(cds_coord)
                print(len(cds_ovlp))

                #if more than 50% of nt is ovl remove
                if (len(cds_ovlp)/len(cds_coord)) > 0.5:
                    print("> 50%")
                    continue
                    removed+=1
                #if total seq < 30 remove
                elif (len(cds_coord) - len(cds_ovlp)) < 30:
                    print("<30")
                    continue
                    removed+=1

                #cut out the ovlp region
                else:
                    org_seq = (ref_cds_handle.readline()).strip("\n")

                    #Open cut_cds fasta
                    cut_cds_path = "/home/sareh/new_cds/cut_cds"
                    cut_cds_file_name = ("cutter_cds_" + ref_cds + ".fa")
                    cut_cds_path = os.path.join(cut_cds_path,cut_cds_file_name)

                    if os.path.isfile(cut_cds_path):
                        handle_cutcds = open(cut_cds_path,"r")

                        ### DO I REMOVE GAPS BEFORE TRIMMING?
                        #nogap_org_seq = org_seq.replace("-","")

                        new_ref_seq = trimed_seq(org_seq,cds_ovlp,start)

                        out_dir = "/home/sareh/ovlp/no_ovlp_cds"
                        out_path = os.path.join(out_dir, ref_cds)
                        out_handle = open(out_path, "w")

                        #write out trimmed ref cds
                        out_handle.write("{}{}\n".format(header,new_ref_seq))

                        #loop through fasta entries
                        queries = convert_fasta(handle_cutcds)
                        for qh, qs in queries:

                            # pairwise alignment
                            qgene_mafft, rgene_mafft = mafft(query=qs,
                            ref=org_seq, trim=False)

                            new_seq = trimed_seq(qgene_mafft,cds_ovlp,start)
                            #write out trimmed cut_cds
                            out_handle.write(">{}\n{}\n".format(qh,new_seq))

                            print(len(qgene_mafft))
                            print(len(rgene_mafft))

    print(removed)
if __name__ == "__main__":
    main()
