# Surfaces

# Table of Contents

1. [Requirements](#requirements)
2. [Usage](#usage)

### Requirements
```
torch >= 1.10.2
tqdm >= 4.62.3
[tmbed](https://github.com/BernhoferM/TMbed) >= 1.0.0
```

### Usage
Reference CDS of genomes can be retrieved from the NCBI Viral Genomes database using the `retrieve_ref_cds_poly.py` script.
```console
python3 polyprotein/retrieve_ref_cds_poly.py --accn fltr_all_ref_accn.csv --dir rerun_outputs/01_ref_cds
```
Nucleotide sequences in fasta files in a directory can be translated to aa sequences using the `nuc_to_aa.py` script
```console
python3 rerun/scripts/nuc_to_aa.py --indir rerun_outputs/01_ref_cds  --outdir rerun_outputs/02_aa_ref_cds
```
After exctracting CDS and translating to amino acids, sequences can be labeled as SE or non-SE using the transmembrane predictor TMbed.
```console
python3 SE/TMbed.py --indir rerun_outputs/02_aa_ref_cds/NC_000858 --outdir rerun_outputs/03_tmbed
```