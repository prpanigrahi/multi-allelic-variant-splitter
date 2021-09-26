# multi-allelic-vcf-splitter
A Python script to process vcf file such that multi-allelic variants are split into individual vcf rows.

# Usage
python3 vcf_split_multi_allelic_variants.py --invcf <input vcf> --outvcf  <output vcf>

# Requirements
- Split INFO and FORMAT fields allele wise.
- For every TAG in INFO and FORMAT, corresponding vcf headers should be present
- Script should be generic enough to handle vcf without format and sample column, single sample vcf file. For multi-sample vcf file consider first sample
- Follow VCF specs to split INFO and FORMAT field
  - IF tag type is 'A': Number of ALT = Number of values
  - IF tag type is 'R': Number of ALT + 1 = Number of values i.e. values for reference allele is also present
  - If tag type is '.': Deduce A or R type based on Number of alleles and Number of values
- Split Genotype as per below
  - If no GT info (., ./., etc), keep as is
  - If GT: 0/0, for all ALT records genotype would be 0/0
  - If GT: 0/1, for first ALT allele, genotype would be 0/1 and for rest it is ./.
  - If GT: 1/0, for first ALT allele, genotype would be 1/0 and for rest it is ./.
  - If GT: 1/1, for first ALT allele, genotype would be 1/2 and for rest it is ./. and so on
  - If GT: 1/2, for first ALT allele, genotype would be 1/. and for second it is ./1 and for rest it is ./.
- Retain phase information (/  or |)



