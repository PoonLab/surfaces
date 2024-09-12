
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

Example, for all genes a virus

```bash
for i in *.fasta; do python3 ../../surfaces/scripts_pipeline/codon_align.py "$i" -o "${i%_step1.fasta}_step2.fasta"; done
```

### Step 3: Manually review the sequence alignment and create the phylogenetic tree

Manually review and remove any problematic sequences from each alignment using Aliview. After removing the sequences, re-align the data and save the alignments with the extension _step3.fasta in the step 3 folder.

Then, run FastTree to build the phylogenetic tree.

For this step, it is necessary to run the FastTree program with all alignments from step 3. Use the following command line
```bash
for f in *.fasta; do  fasttree -nt -quote "$f" >  ${f%.fasta}.nwk ;  done 
```

### Step 4: Prune the phylogenetic trees 

Prune the phylogenetics tree, for this is necessary to run the follow comand  

```bash
for f in *.nwk; do   python3 ../../../surfaces/scripts_pipeline/prunetree.py "$f" > "${f%_step3.nwk}_step4.csv"; done

```
Run the R script named step4_filter.R to generate the graph showing the number of sequences to filter, with the results you get the number of tips.
Then run prunetree.py again and specify the number of tips that were obtained for that protein. For example: 

```bash
python3  ../../../surfaces/scripts_pipeline/prunetree.py ../step3/zika_anchored_capsid_protein_C_step3.nwk  --seq  ../step3/zika_anchored_capsid_protein_C_step3.fasta -t 78 --mode ntips -o  zika_anchored_capsid_protein_C_step4.fasta  --csvfile  zika_anchored_capsid_protein_C_step4.labels.csv
```

