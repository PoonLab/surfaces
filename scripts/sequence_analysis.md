# Measuring selection in virus proteins

Follow these steps in the surfaces Project to consistently analyze our data.

 
## 1. Data collection

### 1.1 Get accession numbers

- To get accession numbers, search your virus from the `taxonomy` browser in NCBI.

- Click your virus species and filter genomes by length (based on the length of your genome). This should give you accession numbers of entries properly classified as your virus species that are complete or close to complete genomes.

- In the right panel, click `Send to`, `Complete Record`. Destination should be `File` and then chose the format `Accession List` and click `Create File`.

 Example of a list with accessions numbers:
```
KU728743.1
KU728742.1
MN125030.1
MN125029.1
```

### 1.2 Download coding sequences

Use `get_all_accns.py` to download the information associated with your genomes.
#### Inputs:
- List with accession numbers
- Email

#### Outputs:
The outputs of this script depends on the type of genome:

**Option 1: Genomes with genes translated independently:**

- `_md.csv`: metadata

- `_aa.fasta`: amino acid sequences

- `_CDS.fasta`: coding sequences (CDSs)

  

Note: We will use the amino acid sequences to cluster the more similar proteins, but the selection analysis will be performed on the CDSs.

  

Example:
``` console
$ python3 get_all_accns.py data/measles_accns.seq user@gmail.com --outfile data/measles
```

**Option 2: Genomes that encode a single polyprotein (run under `--poly` mode):**

- `_md.csv`: metadata

- `_CDSs.fasta`: coding sequences (CDSs)

Note: We will extract the protein encoding sequences by aligning the CDS from the polyprotein to a reference genome in gbk format  with `extract_mat_peptides.py`.

Example:
``` console
$ python3 get_all_accns.py data/potyvitus_accn.seq user@gmail.com --outfile data/potyvirus --poly

```

If you want to obtain the CDS of a virus that encoding by a polyprotein, run the following command in the terminal
Example:
``` console
$ python3  extract_mat_peptides.py  data/dengue_1_genome_reference.gb  data/dengue1_CDS_polyprot.fasta 

```



**Note:** Headers should have the format: accession number, organism name, protein product, strand, and location separated by dashes (`-`) with no spaces.
Header example: `KU728743.1-Measles_virus_genotype_D8-nucleoprotein-1-97_1675`

**Note:**  If you want to obtain the CDS of a virus that encoding by a polyprotein, is necessary to download the archive gbk of genome reference from NCBI  https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/


## 2. Grouping homologous coding sequences

Based on the type of genome there are two different approaches

### 2.1 Grouping independently encoded proteins (3 steps):

Since we can't rely on record annotations to classify our CDSs, we find homologous sequences by calculating k-mer distances between amino acid sequences, and then clustering proteins based on those distances.

### A. Calculate k-mer distances between amino acid sequences

From the multi fasta file with all the amino acid sequences from all CDSs of your virus of interest, calculate a k-mer distance between proteins using `kmer.py`.

**Note:** Even though your search should only contain accession numbers that correspond to full genomes, it could happen that you end up with some partial sequences.

Use the option `--filter` to analyze only those genomes with at least `mean-2` proteins, where mean is the mean number of proteins per genome in your dataset.

#### Inputs:

- Fasta file with amino acid sequences

- Name of your output file

#### Output:

- csv file with kmer distance between sequences

Example:
```console
$ python3 kmer.py measles_feb27_aa.fasta measles_feb27_kmerdist_filtered.csv --filter

```

### B. Extract clusters of homology

Use `hclust.R` to cluster your amino acid sequences based on the k-mer distances.

#### Input:

- csv file with k-mer distance between sequences

In RStudio, provide your working directory, name of file with your k-mer distance, and name of your output file by modify the first few lines:

```R
########################################
# Modify your filenames
########################################
wd <- '/home/laura/Projects/surfaces'
# input file: .csv file with kmer distances, rownames should be fasta headers
kmer.filename <- 'scripts_pipeline/temp_measles/measles_kmerdist.csv'  #input
info.filename <- 'scripts_pipeline/temp_measles/measles_info2.csv'  # output
```

#### Output:

- csv file with cluster classification for each protein and name of the clusters for each protein

**Note**: You can look at successfully clustered proteins and outliers by counting the number of proteins per cluster with `table(info$clusters)`:

```R
>  table(info$clusters)
1  2  3  4  5  6  7  8  9  10
703  699  447  473  698  704  703  703  3  2
```

From this output, you can notice that clusters `9` and `10` only have 3 and 2 proteins assigned, therefore we can consider them outliers.
The number of clusters should be equal to the number of proteins your genome encodes. For mumps, that number is 8. 

### C. Separate sequences in a cluster into independent fasta files
Once you have clustered your data, you can create multiple CDSs with the sequences assigned to the same cluster using `split_CDS_by_cluster.py`
#### Inputs:
- fasta file with all CDSs from all genomes
- `-ci`: csv file with information about proteins assigned to clusters
- `-n`: number of proteins (to decide the maximum number of clusters to analyze --- should be the same as the number of proteins in your genome).
- `-l`: location of your output files

Example
```
python3 split_CDS_by_cluster.py temp_mumps/info/mumps_CDSs.fasta -ci temp_mumps/info/mumps_info2.csv -l temp_mumps/info/
```
#### Outputs:
- (n) number of fasta files of CDSs assigned to the same cluster

### D. (Optional) Visualize the cluster classification on your genome  

Use `plot-genome.py` to check your clustering results.
Each line on the plot correspond to a genome, and the different fragments correspond to genes within the genome.
If the process worked properly, you should see a plot with columns colored consistently across all genomes.

#### Input:
- csv file with proteins assigned to clusters
#### Output:
- pdf with genomes colored based on cluster classification

Example
```console
python3 plot-genome.py measles-protein-clusters-info.csv --outfile measles-genome
```
### 2.2 Grouping coding sequences within a polyprotein:

If your virus encode a single poly protein that is later processed into different proteins, align your CDS against a reference genome with `mat_peptide` annotation with `extract_mat_peptides`:

#### Inputs:

- Genbank file of reference genome with `mat_peptide` annotations

- Fasta file with CDS from the polyprotein

#### Output:

- N fasta files with coding sequences associated with the N proteins in your genome

Example

```console
python3 extract_map_peptides.py data/ref_potyvirus.gb data/potyvirus_CDSs.fasta --label data/poty
```



## 3. Separate proteins and measure selection

The script `selection_by_cluster.py` incorporates the following steps:

- create codon-aware alignments of your CDSs with `mafft`

- construct phylogeny with `fasttree`

- prune your phylogeny to a given length

- estimate selection with FUBAR in `hyphy`

### Inputs:

- fasta files with grouped CDSs

There are four options in this script to consider:

- `--prune`: Prune your trees to reduce their length. For consistency in our analysis we are setting the tree length to `0.5`. The files `.final.mafft.fa` and `.tree` should contain only those sequences that were NOT pruned.

- `--label`: To provide the path and name of your files.

- `--save_before_prune`: For debugging purposes, safe the codon-aware alignment and phylogeny before pruning.

- `--run_sel`: Under this tag the script will measure selection in all clusters of CDSs from your pruned alignments and trees. Note that when this tag is used, the following files will be generated:
### Outputs (per cluster):

1. `.codon_aln.mafft.fasta`, codon-aware aligned classified into the same cluster (after pruning if `--prune` was used)

2. `.cleaned.mafft.fa`, fasta of `mafft` aligned sequences with STOP codons removed by`hyphy`

3. `.tree`, newick file created with `fasttree` (after pruning if `--prune` was used)

4. `FUBAR.json`, json file with selection measurements

5. `.cleaned.mafft.fa.FUBAR.cache`, cache file created during `hyphy` runs

  

### 3.1 Explore your data
Before we run the selection analysis, we need to get familiar with our datasets by exploring the alignments and trees associated to each dataset. 
You can run `selection_by_cluster` using the options `--save_before_prune` and `--save_prune_report` to check whether you can reduce your datasets to a tree with a length around a given target:
```console
python3 selection_by_cluster.py temp_measles/cdss/*.fasta --label temp_mumps/results/measles --prune 1.0 --save_before_prune --save_prune_report
```
The `.report.csv` file should give you an idea of the length of your trees and whether or not they can be pruned to your target length:
|cluster                 |initial_seqs|pruned|inital_tree_lenght|final_seqs|final_tree_length|selection|
|------------------------|------------|------|------------------|----------|-----------------|---------|
|C_protein               |473         |False |0.539             |NA        |NA               |False    |
|fusion_protein          |704         |False |0.796             |NA        |NA               |False    |
|hemagglutinin_protein   |703         |False |0.89              |NA        |NA               |False    |
|large_polymerase_protein|703         |False |0.808             |NA        |NA               |False    |
|matrix_protein          |703         |False |1.094             |NA        |NA               |False    |
|nucleocapsid_protein    |703         |False |0.883             |NA        |NA               |False    |
|phosphoprotein          |699         |False |0.793             |NA        |NA               |False    |
|V_protein               |447         |False |0.726             |NA        |NA               |False    |


### 3.2 Look at your trees and your codon-aware alignments
The former command you will create, for each CDS, a codon aware alignment and a phylogenetic tree:
```
temp_mumps/results$ ls -la *.tree
36 -rw-rw-r-- 1 laura laura 33969 Jul  8 10:22 measles_C_protein.before_prun.tree
56 -rw-rw-r-- 1 laura laura 55303 Jul  8 10:22 measles_fusion_protein.before_prun.tree
60 -rw-rw-r-- 1 laura laura 58134 Jul  8 10:23 measles_hemagglutinin_protein.before_prun.tree
64 -rw-rw-r-- 1 laura laura 63656 Jul  8 10:23 measles_large_polymerase_protein.before_prun.tree
56 -rw-rw-r-- 1 laura laura 54595 Jul  8 10:24 measles_matrix_protein.before_prun.tree
56 -rw-rw-r-- 1 laura laura 56937 Jul  8 10:24 measles_nucleocapsid_protein.before_prun.tree
60 -rw-rw-r-- 1 laura laura 58365 Jul  8 10:24 measles_phosphoprotein.before_prun.tree
36 -rw-rw-r-- 1 laura laura 33342 Jul  8 10:24 measles_V_protein.before_prun.tree
```

```
/temp_measles/results$ ls -la *.codon_aln.mafft.fasta
-rw-rw-r-- 1 laura laura  292902 Jul  8 10:22 measles_C_protein.codon_aln.mafft.fasta
-rw-rw-r-- 1 laura laura 1275172 Jul  8 10:22 measles_fusion_protein.codon_aln.mafft.fasta
-rw-rw-r-- 1 laura laura 1358016 Jul  8 10:22 measles_hemagglutinin_protein.codon_aln.mafft.fasta
-rw-rw-r-- 1 laura laura 4657977 Jul  8 10:23 measles_large_polymerase_protein.codon_aln.mafft.fasta
-rw-rw-r-- 1 laura laura  784046 Jul  8 10:24 measles_matrix_protein.codon_aln.mafft.fasta
-rw-rw-r-- 1 laura laura 1154527 Jul  8 10:24 measles_nucleocapsid_protein.codon_aln.mafft.fasta
-rw-rw-r-- 1 laura laura 1151802 Jul  8 10:24 measles_phosphoprotein.codon_aln.mafft.fasta
-rw-rw-r-- 1 laura laura  463642 Jul  8 10:24 measles_V_protein.codon_aln.mafft.fasta
```

Use `figtree` and `Aliview` or any other programs to visualize the trees identify unusually long branches, and alignments to identify possible frame shifts or other unusual events:
```
temp_measles/results$ figtree measles_V_protein.before_prun.tree 
temp_measles/results$ aliview measles_V_protein.before_prun.tree 
```
### 3.3 Remove potentially problematic sequences from original CDS files
From exploring the trees and codon-aware alignments you will now have an idea of your problematic sequences. 
Use `Aliview` to visualize each one of your CDS files (inputs on step `3.1`) and remove sequences that created long branches on the tree, incomplete sequences, or sequences with a lot of `N` characters in the middle. 

### 3.4 Measure selection
Now that you have manually cleaned your CDS files, you can now proceed to measure selection in all your alignments with the `--run_sel` option on `selection_by_cluster`:

```
python3 selection_by_cluster.py temp_mumps/*.fasta --label temp_mumps/temp/bad --prune 1.0 --save_prune_report --run_sel
```

## 6. Plot selection grid

Use `dnds_grid_plot.R` to visualize your selection estimates for each group of proteins

  

## 7. Measure tree length and store results

In `dnds_grid_plot.R` there is a section to extract the most common name of the proteins in each cluster by loading the `clusters-info.csv` file:

  

```R
md<-read.csv("feb28measles_protein-clusters-info.csv")
g.n<-table(md$gene.name,md$clusters)
prot.name<-c()
for(i  in  1:8){
  m  <-  which.max(g.n[,i])
  prot.name<-append(prot.name,names(m))
}
```

  

Additionally, we can use the following code to calculate the length of the trees from all your proteins (`.tree` created with `selection_by_cluster.py`):

```R
require(ape)
phy.files  <-  Sys.glob("*.tree")
all.trees  <-  lapply(phy.files, function(f) read.tree(f))
t.l  <-  lapply(all.trees, function(t) sum(t$edge.length))
```

Finally, store the information about the proteins in all our clusters in a final table, that we can edit and use as an input for our analysis with all proteins from all viruses:

  

```R
prot.d  <-  data.frame(
  virus =  rep("measles", length(phy.files)),
  cluster =  c(1:length(phy.files)),
  prot.name =  prot.name,
  tree.length =  unlist(t.l),
  is.surface =  rep(NA, length(phy.files))
)
write.csv(prot.d, "measles-surfaces-info.csv", row.names =  FALSE)
```

**Note:** we will need to manually edit the column `is.surface` as `TRUE` or `FALSE` based on the knowledge of the virus.
