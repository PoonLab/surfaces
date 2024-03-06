# Retrieve and cluster proteins from a single virus
Follow these steps in the surfaces Project to consistently analyse our data

## 1. Data collection: complete genomes, amino acid sequences, and CDSs
From a list of accession numbers, use `get_all_accns.py` to create three files:
- `_md.csv`: metadata
- `_aa.fa`: amino acid sequences
- `_CDS.fa`: coding sequences (CDSs)
```
$ python3 get_all_accns.py measles_accns.seq aleatorios.lauram@gmail.com --outfile measles_feb27
```
**Note:** We will use the amino acid sequences to cluster the more similar proteins, but the selection analysis will be performed on the CDSs

## 2. Calculate k-mer distances between amino acid sequences
From the multi fasta file with all the amino acid sequences from all CDSs of your virus of interest, calculate a k-mer distance between proteins using `kmer.py`. 

**Note:** Even though your search should only contain accession numbers that correspond to full genomes, it could happen that you end up with some partial sequences. 
Use the option `--filter` to analyze only those genomes with at least `mean-2` proteins, where mean is the mean number of proteins per genome in your dataset. 

- Input: Fasta file with amino acid sequences
- Output: csv file with kmer distance between sequences

```console
$ python3 kmer.py measles_feb27_aa.fasta measles_feb27_kmerdist_filtered.csv --filter
```

## 3. Extract clusters of homology
Use `hclust.R` to cluster your amino acid sequences based on the kmer distances.

- Input: csv file with kmer distance between sequences
- Output: csv file with proteins assigned to clusters

In RStudio, modify the first few lines to set your working directory and read the kmer distance file
```R
setwd('/home/laura/surfaces/data/measles')
km <- read.csv('measles_distance_filtered.csv', header=F, row.names=1)
```
Additionally, modify the last line to specify the name of your output file:
```R
write.csv(info, 'measles_protein-clusters-info.csv')
```

## 4. Visualize the cluster classification on your genome
Use `plot-genome.py` to check your clustering results.
Each line on the plot correspond to a genome, and the different fragments correspond to genes within the genome. 
If the process worked properly, you should see a plot with columns coloured conssitently across all genomes. 

- Input: csv file sith proteins assigned to clusters
- Output: pdf with genomes colored based on cluster classification

```console
python3 plot-genome.py measles-protein-clusters-info.csv --outfile measles-genome
```

## 5. Separate proteins and measure selection
Use `selection_by_cluster.py` to divide your CDS data into multiple fasta files with coding sequences that correspond to the same protein. 

**Caution:** If you use the the tag `--run_sel` a lot of files with be generated!! 

- Inputs: 
1. csv file sith proteins assigned to clusters
2. Fasta file with all CDSs from all genomes
- Outputs (per cluster):
1. `.fa`, fasta file with CDSs classified into the same cluster
2. `.fa.mafft`, fasta file of `mafft` aligned sequences
3. `.cleaned.fa`, fasta file with STOP codons removed by`hyphy`
4. `.cleaned.tree`, newick file created with `fasttree`
5. `FUBAR.json`, json file with selection measurments 

**Note:** Use the tag `--label` to provide the path and name of your files. Use the tag `--no_prots` to decide the maximum number of clusters to analyze (should be the same as the number of proteins in your genome)
```console
python3 selection_by_cluster.py measles-protein-clusters-info.csv measles_CDSs.fasta --label measles --n_prots 8 --run_sel
```

## 6. Plot selection grid
Use `dnds_grid_plot.R` to visualize your selection estimates for each cluster of proteins
