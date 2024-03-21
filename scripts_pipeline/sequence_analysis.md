# Retrieve and cluster proteins from a single virus
Follow these steps in the surfaces Project to consistently analyse our data. 


## 1. Data collection: complete genomes, amino acid sequences, and CDSs
To get accession numbers, search your virus from the `taxonomy` browser in NCBI. 
Click your virus species and filter genomes by length (based on the length of your genome).
This should give you accession numbers of entries properly classified as your virus species, that are complete genomes.

From the list of accession numbers, use `get_all_accns.py` to create three files:
- `_md.csv`: metadata
- `_aa.fa`: amino acid sequences
- `_CDS.fa`: coding sequences (CDSs)
```
$ python3 get_all_accns.py measles_accns.seq aleatorios.lauram@gmail.com --outfile measles_feb27
```
**Notes:** 
- If your virus genome encodes a polyprotein, you can use the argument `--poly` to get features of `map_peptides` instead of `CDS`
- We will use the amino acid sequences to cluster the more similar proteins, but the selection analysis will be performed on the CDSs

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

## 7. Measure treee length and store results
In `dnds_grid_plot.R` there is a section to extract the most common name of the proteins in each cluster by loading the `clusters-info.csv` file:

```R
md<-read.csv("feb28measles_protein-clusters-info.csv")
g.n<-table(md$gene.name,md$clusters)
prot.name<-c()
for(i in 1:8){
  m <- which.max(g.n[,i])
  prot.name<-append(prot.name,names(m))
}
```

Additionally, we can use the following code to calculate the length of the trees from all your proteins (`.tree` created with `selection_by_cluster.py`):
```R
require(ape)
phy.files <- Sys.glob("*.tree")
all.trees <- lapply(phy.files, function(f) read.tree(f))
t.l <- lapply(all.trees, function(t) sum(t$edge.length))
```
Finally, store the information about the proteins in all our clusters in a final table, that we can edit and use as an input for our analysis with all proteins from all viruses:

```R
prot.d <- data.frame(
              virus = rep("measles", length(phy.files)),
              cluster = c(1:length(phy.files)),
              prot.name = prot.name,
              tree.length = unlist(t.l),
              is.surface = rep(NA, length(phy.files))
                   ) 

write.csv(prot.d, "measles-surfaces-info.csv", row.names = FALSE)
```
**Note:** we will need to manually edit the column `is.surface` as `TRUE` or `FALSE` based on the knowledge of the virus.

# Determining number of sequences and time scale

## 1. Simulate a random tree with 100 tips:
```R
require ape
r.tree <- rtree(100, rooted = TRUE,
                tip.label = NULL)
write.tree(r.tree, file = "random_tree_100.nwk", append = FALSE,
           digits = 10, tree.names = FALSE)

```
## 2. Simulate sequences with INDELible
Use `create_control_indelible.py` to create a `control.txt` file:
```
python3 create_control_indelible.py random_tree_100.nwk -o random_sim -r -s 0.5
```
Note: when using the option `-r`, the script runs indelible from the control file created. 
Run output will be stored at `random_sim`

Alternatively, run indelible from terminal:
```
indelible control.txt
```
## 3. Measure selection with HyPhy
```
hyphy
1
8
output_file.fas  # sequences from indelible
random_tree100.nwk  # tree from ape
```

## Alternativelyu: run 2 and 3 in a single script for multiple tree lengths
To run indelible and measure selection, all at once, use the script `indelible_plus_selection.py`.
This script uses the function `run_selection_pipeline` from `selection_by_cluster.py`.
```
python3 indelible_plus_selection.py ../data_find_tree_length/random_tree_100.nwk -p ../data_find_tree_length/ -r 1
```

## 4. Compare selection measurments with simulation inputs
Simulation inputs for hyphy were:
```
omegas = [
0.0510204, 0.1530612, 0.2551020, 0.3571428, 0.4591836, 0.5612244, 0.6632652, 0.7653060, 0.8673468,
0.9693876, 1.0714284, 1.1734692, 1.2755100, 1.3775508, 1.4795916, 1.5816324, 1.6836732, 1.7857140,
1.8877548, 1.9897956, 2.0918364, 2.1938772, 2.2959180, 2.3979588, 2.4999996, 2.6020404, 2.7040812,
2.8061220, 2.9081628, 3.0102036, 3.1122444, 3.2142852, 3.3163260, 3.4183668, 3.5204076, 3.6224484,
3.7244892, 3.8265300, 3.9285708, 4.0306116, 4.1326524, 4.2346932, 4.3367340, 4.4387748, 4.5408156,
4.6428564, 4.7448972, 4.8469380, 4.9489788, 5.0510196]

prop = [
    1.063764e-01, 1.464867e-01, 1.401630e-01, 1.223917e-01, 1.023007e-01, 8.333079e-02, 6.673162e-02, 
    5.279576e-02, 4.139382e-02, 3.222714e-02, 2.495008e-02, 1.922789e-02, 1.476163e-02, 1.129629e-02,
    8.620596e-03, 6.562952e-03, 4.985996e-03, 3.780965e-03, 2.862470e-03, 2.163923e-03, 1.633688e-03,
    1.231906e-03, 9.279287e-04, 6.982665e-04, 5.249686e-04, 3.943506e-04, 2.960034e-04, 2.220245e-04,
    1.664245e-04, 1.246710e-04, 9.333909e-05, 6.984377e-05, 5.223631e-05, 3.904912e-05, 2.917804e-05,
    2.179306e-05, 1.627075e-05, 1.214322e-05, 9.059520e-06, 6.756626e-06, 5.037501e-06, 3.754636e-06,
    2.797655e-06, 2.084012e-06, 1.551998e-06, 1.155507e-06, 8.600995e-07, 6.400653e-07, 4.762155e-07]
```

If `[printrates] TRUE`, it will print a `RATES.txt` file like:
```
Site	Class	Partition	Inserted?
1	3	1
2	3	1
3	9	1
4	8	1
5	4	1
6	4	1
7	16	1
8	6	1
9	2	1
10	4	1
11	0	1
12	9	1
```

Now I need to map each site to the omega value it corresponds to, and then I can compare site by site against fubar output

**Note:** HyPhy doesn't like INDELible output when `indel rate = 0`, so I re-aligned the file with `mafft` and now it runs properly.

**Note:** To import the `RATES.txt` files to R from the INDELible outpit, we need to remove the first 9 lines with information about the simulation:
```
sed -i '1,9d' *_RATES.txt 
``` 

