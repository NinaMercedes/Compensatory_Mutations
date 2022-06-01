# Compensatory_Mutations
A systematic framework to identify compensatory mutations from _M. tuberculosis_ whole genome sequences. 
## Packages
Implemented in R version 3.6.1. Requirements include BCFtools (v1.9) and R packages: lavaan (v0.6-10), data.table, optparse, and dplyr.

## The Framework
### Stage 1- Pre-process Data (pre_process.R)
This stage takes gene sequence variations in Variant Call Format (vcf.gz), removes lineage-specific SNPs and samples with missing phenotypes to create a genotype file for downstream analysis. There is also the option to use only missense mutations. Note the metadata file must contain columns 'id' and 'lineage', and name of drug (e.g. 'rifampicin'). Example of how to use below. 
```
Rscript pre_process.R --vcf_gz "tb.vcf.gz" --region "rpoB_region.tsv" ----filter_vcf TRUE --out_vcf "rpoB_filtered.vcf.gz" --make_geno TRUE --out_geno "rpoB.geno" --sample_names "sample_names.txt" --metadata "metadata.csv" --drug "rifampicin" --rem_lineage TRUE --out_geno_nolin "rpoB_nolin.geno" --missense TRUE

Rscript pre_process.R --vcf_gz "tb.vcf.gz" --region "rpoC_region.tsv" ----filter_vcf TRUE --out_vcf "rpoC_filtered.vcf.gz" --make_geno TRUE --out_geno "rpoC.geno" --sample_names "sample_names.txt" --metadata "metadata.csv" --drug "rifampicin" --rem_lineage TRUE --out_geno_nolin "rpoC_nolin.geno" --missense TRUE
```
### Stage 2- Analysis (assoc_test.R and med_test.R)
This stage uses genotype files created in stage 1 and performs association and mediation tests to identify putative compensatory mutations. Example of how to use below. 
```
Rscript assoc_test.R --geno_a "rpoB_nolin.geno" --geno_b "rpoC_nolin.geno" --assoc_test TRUE --chi_file "rpoB_rpoC_chi.csv" --metatdata "metadata.csv" --drug "rifampicin"

Rscript med_test.R --geno_a "rpoB_nolin.geno" --geno_b "rpoC_nolin.geno" --med_test TRUE --chi_file "rpoB_rpoC_chi.csv"  --med_file "rpoB_rpoC_med.csv" --metatdata "metadata.csv" --drug "rifampicin"
```
### Stage 3- Interpretation (interpret_results.R)
This stage interprets the results from stage 2 analyses. How decisions are made are shown in summary figure below. The information files for snps in region 1 should contain the headers 'snp1_pos', 'snp1_ref_allele', 'snp1_alt_allele,	'snp1_amino_change', and the equivalent headers for region 2.
 Example of how to use below. 
```
Rscript assoc_test.R --interpret_results TRUE --med_file "rpoB_rpoC_med.csv" --info_file_r1 "rpoB_info.csv" --info_file_r2 "rpoC_info.csv"  --out_file "significant_mediation_rpoB_rpoC.csv" --chi_only "significant_chi_only_rpoB_rpoC.csv"
```
## Overview

![Supplementary File 1](https://user-images.githubusercontent.com/78304174/171453424-388f8928-c579-4090-b5e4-b1a34017da84.jpg)

