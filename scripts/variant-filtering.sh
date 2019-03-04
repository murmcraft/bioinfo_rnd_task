#!/bin/bash

# Variant filtering

# User input variables
INVCF=$1 # Path to input .vcf.gz file
OUT_PREFIX=$2 # Output file prefix
FASTA=$3 # Reference genome .fasta file
OMNI=$4
G1000=$5 
DBSNP=$6
MILLS=$7

# Create an output directory
DIR=variant-filtering 
rm -rf $DIR 
mkdir -p $DIR
cd $DIR

# Hard thresholds fír filtering
# Run GATK default hard thresholds for SNPs
gatk SelectVariants \
    -V ${OUT_PREFIX}_VQSR.vcf.gz \
    -O ${OUT_PREFIX}_SNPs.vcf.gz \
    -R ${FASTA} \
    --select-type-to-include SNP 

gatk VariantFiltration \
    -V ${OUT_PREFIX}_SNPs.vcf.gz \
    -O ${OUT_PREFIX}_SNPs_filters.vcf.gz \
    -R ${FASTA} \
    --filter-name "GATK-defaults-QD" --filter-expression "QD < 2.0" \
    --filter-name "GATK-defaults-MQ" --filter-expression "MQ < 40.0" \
    --filter-name "GATK-defaults-FS" --filter-expression "FS > 60.0" \
    --filter-name "GATK-defaults-SOR" --filter-expression "SOR > 3.0" \
    --filter-name "GATK-defaults-MQRankSum" --filter-expression "MQRankSum < -12.5" \
    --filter-name "GATK-defaults-ReadPosRankSum" --filter-expression "ReadPosRankSum < -8.0" 

# Run GATK default hard threholds for INDELs
gatk SelectVariants \
    -V ${OUT_PREFIX}_VQSR.vcf.gz \
    -O ${OUT_PREFIX}_INDELs.vcf.gz \
    -R ${FASTA} \
    --select-type-to-include INDEL 

gatk VariantFiltration \
    -V ${OUT_PREFIX}_INDELs.vcf.gz \
    -O ${OUT_PREFIX}_INDELs_filters.vcf.gz \
    -R ${FASTA} \
    --filter-name "GATK-defaults-QD" --filter-expression "QD < 2.0" \
    --filter-name "GATK-defaults-MQ" --filter-expression "MQ < 40.0" \
    --filter-name "GATK-defaults-FS" --filter-expression "FS > 200.0" \
    --filter-name "GATK-defaults-SOR" --filter-expression "SOR > 10.0" \
    --filter-name "GATK-defaults-ReadPosRankSum" --filter-expression "ReadPosRankSum < -20.0"

# Merge files together into a single VCF
bcftools concat --allow-overlaps \
    ${OUT_PREFIX}_SNPs_filters.vcf.gz \
    ${OUT_PREFIX}_INDELs_filters.vcf.gz \
    -Oz -o ${OUT_PREFIX}_filters.vcf.gz
bcftools index -t ${OUT_PREFIX}_filters.vcf.gz

# Generate a plottable table from the VCF
gatk VariantsToTable \
    -V ${OUT_PREFIX}_filters.vcf.gz \
    -O ${OUT_PREFIX}_filters.table \
    -F CHROM -F POS -F REF -F ALT \
    -F VQSR \
    -F TYPE -F DP -F MQ -F QD -F SOR \
    -F MQRankSum -F ReadPosRankSum