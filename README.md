
## Pipeline

### 1.1: Identify records from NCBI

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

**Option 1: Genomes with genes translated independently:**

Use `get_all_accns.py` to download the information associated with your genomes.

Inputs:
- List with accession numbers
- Email

Outputs:
The outputs of this script depends on the type of genome:
- `_md.csv`: metadata
- `_aa.fasta`: amino acid sequences (DEPRECATED)
- `_CDS.fasta`: coding sequences (CDSs)

Example:
``` console
$ python3 get_all_accns.py data/measles_accns.seq user@gmail.com --outfile data/measles
```

Next, we need to sort the CDSs for the different gene products (proteins).
Locate a well-annotated reference genome and obtain the Genbank file for that record.
Run the script `sortCDS.py` to partition the CDS file produced by `get_all_accns.py` into multiple FASTA files.

Inputs:
- Genbank file of reference genome
- File containing all CDS records

Output:
- Automatically opens multiple FASTA files to write different CDSs

**Option 2: Genomes encoding a single polyprotein:**

To download sequences from NCBI for genomes with a single open reading frame encoding all proteins, use the `get_all_accns.py` script with the `--poly` flag:
```bash
python3 ../../surfaces/scripts/get_all_accns.py --prefix zika zika.seq --poly hcastelans@gmail.com
```

Next you need to run the `extract_mat_peptides.py` script to extract the CDS.
In this case, it is necessary to use a reference genome in Genbank format to extract the CDS.
```bash
python3 ../../surfaces/scripts/extract_mat_peptides.py '../../surfaces_data/zika/zika/sequence.gb' '/home/hugocastelan/Documents/projects/surfaces_data/zika/zika/zika.seq_CDSs_polyprot.fasta'
```

### Step 2 : Align sequences 

Use the codon_align.py script to align the sequences in codon in frame, example for only one gene. 
```bash
python3 ../../surfaces/scripts/codon_align.py zika_nonstructural_protein_NS4B_step1.fasta -o zika_nonstructural_protein_NS4B_step2.fasta
```

Example, for all genes of a virus
```bash
for i in *.fasta; do python3 ../../surfaces/scripts/codon_align.py "$i" -o "${i%_step1.fasta}_step2.fasta"; done
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
for f in *.nwk; do python3 scripts/prunetree.py "$f" > "${f%_step3.nwk}_step4.csv"; done
```
Then run the R script named `step4_filter.R`  to generate the graph of trend of tree length decay with decreasing number of tips, with the graph you get the number of tips are neccesary to prune. 
NOTE: If length of entire tree is below some threshold (0.5) then abandon alignment (stop here) 

### Step 5: Prune and create phylogenetic trees 
Then run prunetree.py again and specify the number of tips for each protein to necessary to prune accoriding with the results of step 4. For example: 
```bash
python3  ../../surfaces/scripts/prunetree.py ../step3/zika_anchored_capsid_protein_C_step3.nwk  --seq  ../step3/zika_anchored_capsid_protein_C_step3.fasta -t 100 --mode ntips -o  zika_anchored_capsid_protein_C_step5.fasta --csvfile zika_anchored_capsid_protein_C_step5.labels.csv
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
python3 ../../surfaces/scripts/fubar.py  zika_protein_2K_step5.fasta  zika_protein_2K_step5.fubar.csv
```
For batch processing, you can use something like the following:
```bash
for f in data/IBV/*_step5.cds.fa; do python scripts_pipeline/fubar.py $f "${f%.cds.fa}"; done
```
