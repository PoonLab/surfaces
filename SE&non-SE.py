import re
import os
from Bio import SeqIO


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--in', type=str,
                       help='input fasta file of protein seqeunce of each gene')
                       # works with ref_cds (one fasta) or cut_cds (multiple cds)
    parser.add_argument('--out', type=str,
                       help='csv with all the information ')
                       # directoried aa | tmbed | TMHMM
    return parser.parse_args()
  
  
  
  def main()
    args = parse_args()
    
    for file in os.listdir(dir):
        acc <- gsub(".FUBAR.json", "", args.in)
        in_path = os.path.join(dir,accn) 
        aa_path = os.path.join(args.out, "aa", accn)
        
        in_handle = open(in_path,"r")
        out_handle = open(out_path, "w")
        
        # nuc (in) -> aa (out)
        for record in SeqIO.parse(in_handle", "fasta"):
            #print(record.description)
            #print(record.seq)
            record.seq = record.seq.translate() # astricks at end instead of stop codon 
        
            # write out translated fasta file 
            r=SeqIO.write(record, aa_handle, 'fasta')
            if r!=1: print('Error while writing sequence:  ' + seq_record.id)
      
        sleep(2)
       
        out_tmbed = os.path.join(out_dir, "tmbed", accn)
        
        # run tmbed
        run = ("python3 -m tmbed predict -f {} -p {}".format(aa_path ,out_tmbed))
        os.system(run)
        
        # analyse tmbed output 
        
        # accn | number of each alphabet | num aa | 
        

       
if __name__ == '__main__':
    main()  
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
