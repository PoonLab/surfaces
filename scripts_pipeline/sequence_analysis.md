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

Use `get_all_accns.py` to dowload the information associated with your genomes.
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

Note: We will extract the protein encoding sequences by aligning the CDS from the polyprotein to a reference genome with `extract_mat_peptides.py`.

Example:
``` console
$ python3 get_all_accns.py data/potyvitus_accn.seq user@gmail.com --outfile data/potyvirus --poly

```

If you want to obtain the CDs of a virus that codes for a polyprotein, run the following command in the terminal
Example:
``` console
$ python3  extract_mat_peptides.py  data/dengue_1_genome_reference.gb  data/dengue1_CDS_polyprot.fasta 

```



**Note:** Headers should have the format: accession number, organism name, protein product, strand, and location separated by dashes (`-`) with no spaces.
Header example: `KU728743.1-Measles_virus_genotype_D8-nucleoprotein-1-97_1675`


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

In RStudio, modify the first few lines to set your working directory and read the k-mer distance file:

```R

setwd('/home/laura/surfaces/data/measles')

km  <-  read.csv('measles_distance_filtered.csv', header=F, row.names=1)

```

#### Output:

- csv file with cluster classification for each protein

Modify the last line to specify the name of your output file:

```R
write.csv(info, 'measles_protein-clusters-info.csv')
```
**Note**: From line `18`, you should be able to get the number of sequences per cluster. For example:

```R
>  table(clus)
clus
1  2  3  4  5  6  7  8  9  10
703  699  447  473  698  704  703  703  3  2
```

From this output, you can notice that clusters `9` and `10` are outliers.
Thus, when measuring selection you can omit analysis for these clusters by setting the number of proteins (`--n_prots`) to `8`.

### C. (Optional) Visualize the cluster classification on your genome  

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

**From protein clustered with `hclust.R`:**

- csv file with proteins assigned to clusters

- fasta file with all CDSs from all genomes

**From CDSs with already grouped proteins:**

- fasta files with grouped CDSs

There are four options in this script to consider:

- `--prune`: Prune your trees to reduce their length. For consistency in our analysis we are setting the tree length to `0.5`. The files `.final.mafft.fa` and `.tree` should contain only those sequences that were NOT pruned.

- `--label`: To provide the path and name of your files.

- `--no_prots` to decide the maximum number of clusters to analyze (should be the same as the number of proteins in your genome).

- `--save_before`: For debbuding purposes, safe the codon-aware alignment and phylogeny before pruning.

- `--run_sel`: Under this tag the script will measure selection in all clusters of CDSs from your pruned alignments and trees. Note that when this tag is used, the following files will be generated:
### Outputs (per cluster):

1. `.codon_aln.mafft.fasta`, codon-aware aligned classified into the same cluster (after pruning if `--prune` was used)

2. `.cleaned.mafft.fa`, fasta of `mafft` aligned sequences with STOP codons removed by`hyphy`

3. `.tree`, newick file created with `fasttree` (after pruning if `--prune` was used)

4. `FUBAR.json`, json file with selection measurements

5. `.cleaned.mafft.fa.FUBAR.cache`, cache file created during `hyphy` runs

  

### Example if CDSs where clustered with `hclust.R`:

```console
python3 selection_by_cluster.py temp_mumps/mumps_CDSs.fasta -ci temp_mumps/mumps_clusters_info.csv --label temp_mumps/mumps --n_prot 9 --run_sel --prune 0.5 -sb
```

### Example if CDSs where grouped with `extract_mat_peptide.py`:
```console
python3 selection_by_cluster.py temp_hepa/*.fasta --label temp_hepa/ --prune 0.5 -sb --run_sel
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
