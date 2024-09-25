
## Pipeline

### Step 1: Download sequences from NCBI

#### Option 1 

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

#### 1.2 Download coding sequences

Use `get_all_accns.py` to download the information associated with your genomes.
#### Inputs:
- List with accession numbers
- Email

##### Outputs:
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

2. Grouping homologous coding sequences

Based on the type of genome there are two different approaches

### 1.2 Grouping independently encoded proteins (3 steps):

Since we can't rely on record annotations to classify our CDSs, we find homologous sequences by calculating k-mer distances between amino acid sequences, and then clustering proteins based on those distances.

#### A. Calculate k-mer distances between amino acid sequences

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

#### B. Extract clusters of homology

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

#### C. Separate sequences in a cluster into independent fasta files
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
##### Outputs:
- (n) number of fasta files of CDSs assigned to the same cluster

#### D. (Optional) Visualize the cluster classification on your genome  

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

#### Option 2 

To download sequences from NCBI, use the `get_all_accns.py` script:

```bash
python3 ../../surfaces/scripts_pipeline/get_all_accns.py --prefix zika zika.seq --poly hcastelans@gmail.com
```

If your virus is encoded by a polyprotein, you need to run the extract_mat_peptides.py script to extract the CDS. In this case, it is necessary to use a reference genome in GBK format to extract the CDS.

```bash
python3 ../../surfaces/scripts_pipeline/extract_mat_peptides.py '../../surfaces_data/zika/zika/sequence.gb' '/home/hugocastelan/Documents/projects/surfaces_data/zika/zika/zika.seq_CDSs_polyprot.fasta'
```
 
### Step 2 : Align sequences 

Use the codon_align.py script to align the sequences in codon in frame, example for only one gene. 
```bash
python3 ../../surfaces/scripts_pipeline/codon_align.py zika_nonstructural_protein_NS4B_step1.fasta -o zika_nonstructural_protein_NS4B_step2.fasta
```

Example, for all genes of a virus

```bash
for i in *.fasta; do python3 ../../surfaces/scripts_pipeline/codon_align.py "$i" -o "${i%_step1.fasta}_step2.fasta"; done
```

### Step 3: Manually review the sequence alignment and create the phylogenetic tree

Manually review and remove any problematic sequences from each alignment using Aliview. After removing the sequences, re-align the data and save the alignments with the extension _step3.fasta in the step 3 folder.

Then, run FastTree to build the phylogenetic tree.

For this step, it is necessary to run the FastTree program with all alignments from step 3. Use the following command line
```bash
for f in *.fasta; do fasttree -nt -quote "$f" >  ${f%.fasta}.nwk ;  done 
```

### Step 4: Prune the phylogenetic trees 

Prune the phylogenetics tree, for this is necessary to run the following command:
```bash
for f in *.nwk; do python3 scripts_pipeline/prunetree.py "$f" > "${f%_step3.nwk}_step4.csv"; done
```
Then run the R script named `step4_filter.R`  to generate the graph of trend of tree length decay with decreasing number of tips, with the graph you get the number of tips are neccesary to prune. 
NOTE: If length of entire tree is below some threshold (0.5) then abandon alignment (stop here) 

### Step 5: Prune and create phylogenetic trees 
Then run prunetree.py again and specify the number of tips for each protein to necessary to prune accoriding with the results of step 4. For example: 
```bash
python3  ../../../surfaces/scripts_pipeline/prunetree.py ../step3/zika_anchored_capsid_protein_C_step3.nwk  --seq  ../step3/zika_anchored_capsid_protein_C_step3.fasta -t 78 --mode ntips -o  zika_anchored_capsid_protein_C_step5.fasta --csvfile zika_anchored_capsid_protein_C_step5.labels.csv
```

For batch processing, you can run something like the following:
```bash
for f in data/HCV/HCV1a_*_step3.nwk; do python scripts_pipeline/prunetree.py -t 100 \
--mode ntips --seq "${f%.nwk}.fasta" -o "${f%_step3.nwk}_step5.fasta" --csvfile \
"${f%_step3.nwk}_step5.labels.csv" $f; done
```

### Step 6: Selection analysis 

Run the script `fubar.py` to run the selection analysis and then plot the fingerprints with the script `fingerprint_dnds_plot.R` 
```bash
python3 ../../../surfaces/scripts_pipeline/fubar.py  zika_protein_2K_step5.fasta  zika_protein_2K_step5.fubar.csv
```
For batch processing, you can use something like the following:
```bash
for f in data/IBV/*_step5.cds.fa; do python scripts_pipeline/fubar.py $f "${f%.cds.fa}"; done
```
