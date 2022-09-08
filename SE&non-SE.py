import re
import os

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--in', type=str,
                       help='input fasta file of protein seqeunce of each gene')
    parser.add_argument('--out', type=str,
                       help='csv with all the information ')
    return parser.parse_args()
  
  
  
  def main()
    args = parse_args()
    
    for file in os.listdir(dir):
        acc <- gsub(".FUBAR.json", "", args.in)
        in_path = os.path.join(dir,accn) 
        out_path = os.path.join(args.out,accn)
        
        in_handle = open(in_path,"r")
        out_handle = open(out_path, "w")
        
        # nuc (in) -> aa (out)
        
        
        
        sleep(2)
       
        
        # run tmbed
        run = ("python3 -m tmbed predict -f {} -p {}".format(in_path ,out_path))
        os.system(run)
        
        # analyse tmbed output 
        
        # accn | number of each alphabet | num aa | 
        

       
if __name__ == '__main__':
    main()  
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
