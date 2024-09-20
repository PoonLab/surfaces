
## Surfaces Pipeline

### Step 1: Download sequences from NCBI

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

### Step 6: Selection analysis 

Run the script `fubar.py` to run the selection analysis and then plot the fingerprints with the script `fingerprint_dnds_plot.R` 
```bash
python3 ../../../surfaces/scripts_pipeline/fubar.py  zika_protein_2K_step5.fasta  zika_protein_2K_step5.fubar.csv
```
