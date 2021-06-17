#cutter iterator mpi
#python3 cutter.py /home/sareh/surfaces/retrieve_sequences/cutter_test/cutter_ref_example.txt
# /home/sareh/data/Pruned_sequence/Pruned_nuc_NC_001547
# /home/sareh/surfaces/retrieve_sequences/cutter_test/output_NC_00154 --echo
# > /home/sareh/surfaces/retrieve_sequences/cutter_test/NC_001547_score.txt


from mpi4py import MPI
import argparse
from glob import glob
import subprocess
import os

the_world = MPI.COMM_WORLD
my_number = the_world.Get_rank()
total_number = the_world.Get_size()

#parser = argparse.ArgumentParser(description ='path in BEVi to /data/nuc_refseq_CDS')
#parser.add_argument('directory', type =str, help='Path to Directory of virus specied directories with fasta files for each CDS')
#args = parser.parse_args()

directory = '/home/sareh/cutter/rev_cutter/rev_nuc_refseq_CDS'
count = 0
for dir in os.listdir(directory):
    if 'NC_' in dir:
        #path = ("{}/{}".format(directory,dir)) #print(dir) #NC_044047
        path = os.path.join(directory, dir)
        for file in os.listdir(path):
            count += 1 
            if count % total_number == my_number:
                #print(file) #NC_044047_YP_009679042.1
                ref = ("{}/{}/{}".format(directory, dir, file))
                fasta = ("/home/sareh/data/Pruned_sequence/Pruned_nuc_{}".format(dir))
                outfile = ("/home/sareh/cutter/rev_cutter/cutter_cds/rev_cutter_{}".format(file))
                score = ("/home/sareh/cutter/rev_cutter/rev_cutter_scores/{}".format(file))

            	#f.write("python3 cutter.py {} {} {} --echo > {}\n".format(ref,fasta,outfile,score))
                #result=subprocess.check_call(["python3","cutter.py",ref,fasta,outfile,"--echo"])
                with open(score, 'w') as handle:
                    result = subprocess.run(["python3", "cutter.py", ref, fasta, outfile, "--echo"], stdout=subprocess.PIPE)
                    handle.write(result.stdout.decode())
                break

