# Evolutionary fingerprinting of virus proteins: Scripts

## Data collection
* `codon_align.py` - translate a codon alignment into amino acid sequences, perform a multiple sequence alignment (calling out to MAFFT), and then apply the resulting alignment back to the original nucleotide sequences
* `concat_genes.py` - used to deal with spliced exons in IAV and IBV data, for example M2
* `extract_mat_peptides.py` - use Genbank annotations from a reference genome to partition a polypeptide into individual gene products
* `get_all_accns.py` - given a list of Genbank accession numbers, retrieve full genome sequences, coding sequences and metadata using the Entrez API
* `sortCDS.py` - use a reference Genbank record to classify CDS records by pairwise alignment
* `split_genes.py` - split a FASTA file containing two sets of non-homologous gene sequences.  Used for IAV, where one genome segment can encode two genes, e.g., M1 and M2

## Sensitivity analysis
* `create_control_indelible.py` - generate a control file for running a simulation in [INDELible](https://academic.oup.com/mbe/article/26/8/1879/980884) based on a template - used for analysis of tree length requirements
* `indelible_plus_selection.py` - use INDELible to simulate sequences under different tree lengths to evaluate its impact on measuring selection.  Requires `create_control_indelible.py`
* `indelible_vs_fubar.R` - Process and visualize results from `indelible_plus_selection.py`

## Selection analysis
* `dnds.py` - a Python wrapper for the single likelihood ancestor counting (SLAC) method in [HyPhy](https://www.hyphy.org)
* `fubar.py` - Python wrapper for fast unconstrained Bayesian approximation (FUBAR) method in HyPhy
* `mpi_fubar.py` - an MPI-enabled driver script for `fubar.py`
* `plot-dnds.R` - visualize gene-wide dN/dS estimates - runs some generalized linear model analyses to compare between groups

## Downsampling
* `prunetree.py` - progressively remove the shortest terminal branches in a tree until we reach a target tree length
* `sample_codons.py` - sample codon sites from an alignment at random with or without (default) replacement.
* `step4_batch.py` - batch process a set of Newick tree files with `prunetree.py` (automates input and output file naming)
* `step4_filter.R` - visualize the decay of tree length with the progressive removal of the shortest terminal branches.  Locates an "elbow" if user wants to choose a sample size by that criterion (deprecated).
* `subtree_picker.py` - Downsample sequences in a FASTA file by progressively growing the subtree containing a terminal branch specified by the user (by label) until a target total length is obtained.

## Fingerprinting
* `earthmover.R` - calculates the Wasserstein distance matrix for evolutionary fingerprints extracted using 
* `get_fingerprints.R` - extract the evolutionary fingerprint from the FUBAR JSON output and save the entire set in an RData file - also generate heatmap visualizations
* `mds.R` - use multidimensional scaling to visualize Wasserstein distance matrices computed in `earthmover.R`.  Calculate centroids for replicate data sets.  Attempts some SVM classification.
* `wkkn.R` - R script implementing weighted k-nearest neighbours - deprecated by another version.

## Miscellaneous
* `attach_coords.py` - generate a consensus sequence for each protein-coding gene alignment, do pairwise alignment of the consensus to the corresponding reference genome, and extract reference coordinates
* `consensus.py` - generate a consensus sequence from an alignment
* `get_align_stats.py` - fix sequence headers for HyPhy and calculate alignment statistics, including tree length
* `remove_sequences.py` - remove sequences from one FASTA file that are not found in a second FASTA file (pretty trivial)
