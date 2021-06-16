import os
#import minimap2

Thread=1
min=100

infile_directory='/Users/sarehchimeh/Data/Sequences'
refseq_directory='/Users/sarehchimeh/Data/refseq_genome'
file='/Users/sarehchimeh/Data/minimap_iterator.txt'
with open (file,'w+') as f:
    for filename in os.listdir(infile_directory):
        if 'NC_' in filename:
            accn=filename[5:-4]
            outfile = 'minimap_alignment/alignment_{}'.format(accn)
            infile = 'sequences/{}'.format(filename)
            for refname in os.listdir(refseq_directory):
                if accn in refname:
                    refseq ='refseq_genome/{}'.format(refname)
                    #print(infile)
                    #print(refseq)
                    #print(outfile)
                    #print(Thread)
                    #print(min)
                    f.write('python3 minimap2.py' + ' -o ' + outfile + ' --minlen ' + str(min) + ' --ref ' + refseq +' '+ infile+' -f -a\n')
                    #plug in (infile, refseq, outfile, Thread, min) into minimap script here


