**Downsizing**
nuc_genome -> pruned_genome(100)

	   genome file [/home/sareh/2020/sequences/nuc_genome]
           (only ncbi viruses need pruning)

        a. nuc_genome -> msa_genome
           $python3 minimap2.py -o minimap_alignment.fa --ref NC_001802.fa NC_001802 -f -a
           $for f in NC_*; do python3 minimap2.py -o /home/sareh/2020/downsizing/ncbi_genome_msa/$f --ref /home/sareh/2020/sequences/ncbi_ref_nuc_genome/$f /home/sareh/2020/sequences/ncbi_nuc_genome/$f -f -a; done
           *Check number of sequences that aligned (Less)

        b. msa_genome -> fasttree_genome
           double percision FastTree
           $fasttree -nt -gtr -gamma alignment_NC_001802 > fasttree_NC_001802.tre
           $for f in *; do fasttree -nt -gtr -gamma $f > /home/sareh/2020/downsizing/ncbi_genome_fasttree/$f; done

        c. fasttree_genome ---prunetree.py----> pruned_fasttree
           script: fixed_prunetree.py [newick_file target]
           python3 fixed_prunetree.py FastTree_NC_001802.tre 100 > NC_01802_pruned.tre
           python3 prune_fasta.py --tree --fasta --outfile

        d. pruned_fasttree -> pruned_genome
           script: prune_fasta.py [--tree --fasta --outfile]
