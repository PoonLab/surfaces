# Surfaces
### Retrieval of Reference CDS
Reference CDS of genomes can be retrieved from the NCBI Viral Genomes database using the `retrieve_ref_cds_poly.py` script.
```console
python3 polyprotein/retrieve_ref_cds_poly.py --accn fltr_all_ref_accn.csv --dir rerun_outputs/01_ref_cds
```
Nucleotide sequences in fasta files in a directory can be translated to aa sequences using the `nuc_to_aa.py` script
```console
python3 rerun/scripts/nuc_to_aa.py --indir rerun_outputs/01_ref_cds  --outdir rerun_outputs/02_aa_ref_cds
```
