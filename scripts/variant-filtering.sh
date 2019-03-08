#!/bin/bash

# Variant filtering

# User input variables
INVCF=$1 # Path to input .vcf.gz file
OUT_PREFIX=$2 # Output file prefix
FASTA=$3 # Reference genome .fasta file
DBSNP=$4 # Known variants from dbSNP

# Create an output directory
DIR=variant-filtering 
rm -rf $DIR 
mkdir -p $DIR
cd $DIR


# Hard thresholds for filtering
# Run GATK default hard thresholds for SNPs
gatk SelectVariants \
    -V ${INVCF} \
    -O ${OUT_PREFIX}_SNPs.vcf.gz \
    -R ${FASTA} \
    --select-type-to-include SNP 

gatk VariantFiltration \
    -V ${OUT_PREFIX}_SNPs.vcf.gz \
    -O ${OUT_PREFIX}_SNPs_gatkfilters.vcf.gz \
    -R ${FASTA} \
    --filter-name "DP" --filter-expression "DP < 10" \
    --filter-name "QUAL" --filter-expression "QUAL < 50.0" \
    --filter-name "QD" --filter-expression "QD < 5.0" \
    --filter-name "MQ" --filter-expression "MQ < 50.0" \
    --filter-name "FS" --filter-expression "FS > 15.0" \
    --filter-name "SOR" --filter-expression "SOR > 3.0" \
    --filter-name "MQRankSum" --filter-expression "MQRankSum < -2.0" \
    --filter-name "ReadPosRankSum" --filter-expression "ReadPosRankSum < -2.0" 

# Add filtering based on genotype quality
bcftools filter -e 'FMT/GQ < 90' --soft-filter GQ \
    ${OUT_PREFIX}_SNPs_gatkfilters.vcf.gz \
    -Oz -o ${OUT_PREFIX}_SNPs_filters.vcf.gz
bcftools index -t ${OUT_PREFIX}_SNPs_filters.vcf.gz

# Run GATK default hard threholds for INDELs
gatk SelectVariants \
    -V ${INVCF} \
    -O ${OUT_PREFIX}_INDELs.vcf.gz \
    -R ${FASTA} \
    --select-type-to-include INDEL 

gatk VariantFiltration \
    -V ${OUT_PREFIX}_INDELs.vcf.gz \
    -O ${OUT_PREFIX}_INDELs_gatkfilters.vcf.gz \
    -R ${FASTA} \
    --filter-name "DP" --filter-expression "DP < 10" \
    --filter-name "QUAL" --filter-expression "QUAL < 50.0" \
    --filter-name "GATK-defaults-QD" --filter-expression "QD < 5.0" \
    --filter-name "GATK-defaults-MQ" --filter-expression "MQ < 50.0" \
    --filter-name "GATK-defaults-FS" --filter-expression "FS > 10.0" \
    --filter-name "GATK-defaults-SOR" --filter-expression "SOR > 3.0" \
    --filter-name "MQRankSum" --filter-expression "MQRankSum < -2.0" \
    --filter-name "GATK-defaults-ReadPosRankSum" --filter-expression "ReadPosRankSum < -2.0"

# Add filtering based on genotype quality
bcftools filter -e 'FMT/GQ < 90' --soft-filter GQ \
    ${OUT_PREFIX}_INDELs_gatkfilters.vcf.gz \
    -Oz -o ${OUT_PREFIX}_INDELs_filters.vcf.gz
bcftools index -t ${OUT_PREFIX}_INDELs_filters.vcf.gz

# Merge files together into a single VCF
bcftools concat --allow-overlaps \
    ${OUT_PREFIX}_SNPs_filters.vcf.gz \
    ${OUT_PREFIX}_INDELs_filters.vcf.gz \
    -Oz -o ${OUT_PREFIX}_filters.vcf.gz
bcftools index -t ${OUT_PREFIX}_filters.vcf.gz

# Collect metrics for filtered variants
gatk CollectVariantCallingMetrics \
    -I ${OUT_PREFIX}_filters.vcf.gz \
    -O ${OUT_PREFIX} \
    -R ${FASTA} \
    --DBSNP ${DBSNP} 

# Generate a plottable table from the VCF
gatk VariantsToTable \
    -V ${INVCF} \
    -O ${OUT_PREFIX}.table \
    -F CHROM -F POS -F REF -F ALT \
    -F QUAL -F FILTER -F TYPE -F HET \
    -F EVENTLENGTH -F MULTI-ALLELIC -F TRANSITION \
    -F DP -F QD -F FS -F MQ -F SOR \
    -F BaseQRankSum -F MQRankSum -F ReadPosRankSum \
    -GF GQ \
    --show-filtered true \
    --split-multi-allelic true
