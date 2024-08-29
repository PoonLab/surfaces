
## Surface Pipeline

### Step 1: Download Sequences from NCBI

To download sequences from NCBI, use the `get_all_accns.py` script:

```bash
python3 ../../surfaces/scripts_pipeline/get_all_accns.py --prefix zika zika.seq --poly hcastelans@gmail.com
```

If your virus is encoded by a polyprotein, you need to run the extract_mat_peptides.py script to extract the CDS.In this case is necessary to use reference genome to in gbk format to extract the CDS 

```bash
python3 ../../surfaces/scripts_pipeline/extract_mat_peptides.py '/home/hugocastelan/Documents/projects/surfaces_data/zika/zika/sequence.gb' '/home/hugocastelan/Documents/projects/surfaces_data/zika/zika/zika.seq_CDSs_polyprot.fasta'
```

### Step 2 : Align Sequences 

Use the codon_align.py script to align the sequences in codon frame:

```bash
python3 /home/hugocastelan/Documents/projects/surfaces/scripts_pipeline/codon_align.py zika_nonstructural_protein_NS4B_step1.fasta -o zika_nonstructural_protein_NS4B_step2.fasta
```

### Step 3: Manual Sequence Review

Manually review and remove any problematic sequences from each alignment.
