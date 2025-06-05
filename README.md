
## Pipeline

### Step 1: Data collection 
#### 1.1. Identify records from NCBI

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

#### 1.2. Download coding sequences

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
  
Example:
``` console
$ python3 ../../surfaces/scripts/sortCDS.py  rsv_sequence.gb  rsv_sequence.seq_CDSs.fasta
```

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
python3 ../../surfaces/scripts/extract_mat_peptides.py  --label zika 'sequence_zika.gb' '/home/hugocastelan/zika.seq_CDSs_polyprot.fasta' 
```

### Step 2: Sequence alignment

Use the `codon_align.py` script to align your protein encoding sequence in a codon-aware manner:

```bash
python3 ../../surfaces/scripts/codon_align.py zika_nonstructural_protein_NS4B_step1.fasta -o zika_nonstructural_protein_NS4B_step2.fasta
```

- For multiple coding sequences:
```bash
for i in *.fasta; do python3 ../../surfaces/scripts/codon_align.py "$i" -o "${i%_step1.fasta}_step2.fasta"; done
```

### Step 3: Revise sequence alignments and create phylogenies

- Use `Aliview` to visualize your alignments. If required, remove problematic sequences or nucleotides that introduce frameshifts. After removing the sequences, re-align the data and save the alignments with the extension `_step3.fasta` in the step 3 folder.

- Use `FastTree` to build the phylogenetic tree on the manually curated alignments:
```bash
for f in *.fasta; do fasttree -nt -quote "$f" >  ${f%.fasta}.nwk ;  done 
```

### Step 4: Examine decay of tree length and prune tree to a target number of tips

- Use `prunetree.py` to record the effect on tree sizes when progressively down-sampling the sequence alignment by removing the shortest tips of the tree.
```bash
for f in *_step3.nwk; do python3 scripts/prunetree.py "$f" > "${f%_step3.nwk}_step4.csv"; done
```
- Use `step4_filter.R` to create a tree-length-decay plot. Take note of the number of tips to prune per alignment.
- It is necessary to include the "path in quotes"
  
```bash
Rscript step4_filter.R  "/home/hugocastelan/Documents/projects/surfaces_data/dengue/step4/*.csv"  "/home/hugocastelan/Documents/projects/surfaces_data/dengue/step4/figure.png" > step4_filter.csv
```
**NOTE: If length of entire tree is below some threshold (0.5) then abandon alignment (stop here)**

- Use `prunetree.py` one more time, but now provide a target number of tips (`--target` option under `--mode ntips`) for each alignment based on your previous output.  
```bash
python3 prunetree.py measles_C_protein_step3.nwk --seq measles_C_protein_step2.fasta --mode ntips --target 97 --csvfile measles_C_protein_step5.labels.csv --outfile measles_C_protein_step5.fasta
```
For batch processing, use:
```bash
python3 "step4_batch.py" "step4_filter.csv" "*_step3.nwk" 
```
**NOTES:**
* `step4_filter.csv` is a CSV from `step4_filter.R` output
* Script assumes Newick files are in the same directory as corresponding `\*\_step3.fasta` files

### Step 5 (Optional): Prune to target tree length of `1.0` if trees exceed this size
- Case 1: prune down to tree length of `1.0` if doing so does not removes too many tips (change mode to `--mode treelen`). 

- Case 2: prune to `100` tips when you cannot reach a tree length of `1.0` without removing too many sequences.
```bash
for f in data/HCV/HCV1a_*_step3.nwk; do python scripts_pipeline/prunetree.py -t 100 \
--mode ntips --seq "${f%.nwk}.fasta" -o "${f%_step3.nwk}_step5.fasta" --csvfile \
"${f%_step3.nwk}_step5.labels.csv" $f; done
```

### Step 6: Selection analysis 

- Use `fubar.py` to run FUBAR on your pruned alignments and store selection pressures (dN and dS per codon) in a `csv` file.
For example:
```bash
python3 fubar.py  zika_protein_2K_step5.fasta  zika_protein_2K_step5.fubar.csv
```
For batch processing, use:
```bash
for f in *step5.fasta; do python3 scripts/fubar.py "$f" "${f%.fasta}"; done
```

- Create and visualize selection fingerprints with `fingerprint_dnds_plot.R`. 
