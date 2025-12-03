# Evolutionary fingerprinting of virus proteins: Scripts

* `attach_coords.py` - generate a consensus sequence for each protein-coding gene alignment, do pairwise alignment of the consensus to the corresponding reference genome, and extract reference coordinates
* `codon_align.py` - translate a codon alignment into amino acid sequences, perform a multiple sequence alignment (calling out to MAFFT), and then apply the resulting alignment back to the original nucleotide sequences
* `concat_genes.py` - used to deal with spliced exons in IAV and IBV data, for example M2
* `consensus.py` - generate a consensus sequence from an alignment
* `create_control_indelible.py` - generate a control file for running a simulation in [INDELible](https://academic.oup.com/mbe/article/26/8/1879/980884) based on a template - used for analysis of tree length requirements
* `dnds.py` - a Python wrapper for the single likelihood ancestor counting (SLAC) method in [HyPhy](https://www.hyphy.org)
* `earthmover.R` - calculates the Wasserstein distance matrix for evolutionary fingerprints extracted using 
* `extract_mat_peptides.py` - use Genbank annotations from a reference genome to partition a polypeptide into individual gene products
