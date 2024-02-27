# Retrieve and cluster proteins from a single virus
Follow these steps in the surfaces Project to consistently analyse our data

## 1. Data collection: complete genomes, amino acid sequences, and CDSs
From a list of accession numbers, use `get_all_accns.py` to create three files:
- metadata
- amino acid sequences
- coding sequences (CDSs)
```
$ python3 get_all_accns.py measles_accns.seq aleatorios.lauram@gmail.com --outfile measles_feb27
```
**Note:** We will use the amino acid sequences to cluster the more similar proteins, but the selection analysis will be performed on the CDSs

## 2. Calculate k-mer distances between amino acid sequences
From the multi fasta file with all the amino acid sequences from all CDSs of your virus of interest, calculate a k-mer distance between proteins using `kmer.py`. 

**Note:** Even though your search should only contain accession numbers that correspond to full genomes, it could happen that you end up with some partial sequences. Use the option `--filter` to analyse only those genomes with at least `mean - 2` proteins, where mean is the mean number of proteins per genome in your dataset. 

- Input: Fasta file with amino acid sequences
- Output: csv file with kmer distance between sequences

```console
$ python3 kmer.py measles_feb27_aa.fasta measles_feb27_kmerdist_filtered.csv --filter
```

## 3. Extract clusters of homology
Use `hclust.R` to cluster your amino acid sequences based on the kmer distances.

- Input: csv file with kmer distance between sequences
- Output: csv file with proteins assigned to clusters

In RStudio, modify the first few lines to set your working directory, input the kmer distance output, and define the number of proteins in your genome. 
```R
setwd('/home/laura/Projects/ProtClust/scripts')
km <- read.csv('measles_distance_filtered.csv', header=F, row.names=1)
n.prots <- 8 # number of proteins per genome
km <- 1-as.matrix(km)
```
Additionally, modify the last line to specify the name of your output file:
```R
write.csv(info, 'measles_protein-clusters-info.csv')
```

## 4. Visualize the cluster classification on your genome
- Input: csv file sith proteins assigned to clusters
- Output: pdf with genomes colored based on cluster classification

```console
python3 plot-genome.py measles-protein-clusters-info.csv --outfile measles-genome
```