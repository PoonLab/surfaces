
## Surfaces Pipeline

### Step 1: Download Sequences from NCBI

To download sequences from NCBI, use the `get_all_accns.py` script:

```bash
python3 ../../surfaces/scripts_pipeline/get_all_accns.py --prefix zika zika.seq --poly hcastelans@gmail.com
```

If your virus is encoded by a polyprotein, you need to run the extract_mat_peptides.py script to extract the CDS. In this case, it is necessary to use a reference genome in GBK format to extract the CDS.

```bash
python3 ../../surfaces/scripts_pipeline/extract_mat_peptides.py '../../surfaces_data/zika/zika/sequence.gb' '/home/hugocastelan/Documents/projects/surfaces_data/zika/zika/zika.seq_CDSs_polyprot.fasta'
```
 
### Step 2 : Align Sequences 

Use the codon_align.py script to align the sequences in codon in frame, example for only one gene. 
```bash
python3 ../../surfaces/scripts_pipeline/codon_align.py zika_nonstructural_protein_NS4B_step1.fasta -o zika_nonstructural_protein_NS4B_step2.fasta
```

Example, for all genes a virus

```bash
for i in *.fasta; do python3 ../../surfaces/scripts_pipeline/codon_align.py "$i" -o "${i%_step1.fasta}_step2.fasta"; done
```

### Step 3: Manual Sequence Review

Manually review and remove any problematic sequences from each alignment using aliview


### Step 4: Create and prune the phylogenetic trees 

For this step, it is necessary to run the FastTree program with all alignments from step 3. Use the following command line

```bash
for f in *.fasta; do o=$(echo $f | sed 's/.fasta/.nwk/g');   fasttree   -nt -quote  $i > $o ; done
```

or
```bash
for f in *.fasta; do  fasttree -nt -quote "$f" >  ${f%.fasta}.nwk ;  done 
```
Once the phylogenetic trees are built, it is necessary to visualize them in FigTree to check for outliers. If an outlier is found, it should be removed from the alignment in step 3
