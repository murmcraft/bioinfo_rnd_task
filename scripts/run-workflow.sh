#!/bin/bash

# Complete pipeline to run Steps 1. - 5.

# Outputs all variant calls in a VCF file and 
# HTML reports of BAM and variant quality.
# See README.md for details.

# Requirements:
# - BAM file passes validate-sam-file.sh
# - Required reference files at ${WKD}/reference-data/

# User input parameters
WKD=$1 # Path to workdir where input data is and output will be written, e.g. /home
SCRIPTS=$2 # Path to scripts dir, e.g. /scripts
DATASET=$3 # Dataset name, e.g. NA12878
FASTA=$4 # Path to reference genome .fasta file
DBSNP=$5 # Path to known variants in dbSNP as a VCF file
INDELS=$6 # Path to known indels as a VCF file
INTERVALS=$7 # Path to list of non-overlapping genomic intervals for parallelization
SAMPLE_SIZE=$8 # Fraction of data to use for subsetting 

set -e

# Confirm prerequisites
if [[ ! -d "${WKD}" ]]; then
    echo "Error: workdir missing"
    exit 1
fi
if [[ ! -d "${SCRIPTS}" ]]; then
    echo "Error: scripts dir missing"
    exit 1
fi
if [[ ! -s "${FASTA}" ]]; then
    echo "Error: reference genome .fasta missing"
    exit 1
fi
if [[ ! -s "${DBSNP}" ]]; then
    echo "Error: known variant sites VCF missing"
    exit 1
fi
if [[ ! -s "${INDELS}" ]]; then
    echo "Error: known indel sites VCF missing"
    exit 1
fi
if [[ ! -s "${INTERVALS}" ]]; then
    echo "Error: genomic intervals file missing"
    exit 1
fi
if [[ ! ("${SAMPLE_SIZE}" > 0 && "${SAMPLE_SIZE}" < 1) ]]; then
    echo "Error: subsetting fraction has to be above 0 and below 1"
    exit 1
fi

echo "Prerequisites met - processing started"

# Create a dir for run logs
mkdir -p ${WKD}/run-logs

# Step 1 - Mark duplicates
echo "Step 1 - Mark duplicates - Started"
mark-duplicates.sh \
    ${WKD}/${DATASET}.bam \
    ${DATASET}_markdups \
    2>&1 > ${WKD}/run-logs/mark-duplicates.log
echo "Step 1 - Mark duplicates - Finished"
echo "Step 1 outputs stored in: "${WKD}/mark-duplicates

# Step 2 - Base Quality Score Recalibration
echo "Step 2 - Base Quality Score Recalibration - Started"
base-quality-score-recalibration.sh \
    ${WKD}/mark-duplicates/${DATASET}_markdups.bam \
    ${DATASET}_bqsr \
    ${FASTA} \
    ${DBSNP} \
    ${INDELS} \
    2>&1 > ${WKD}/run-logs/base-quality-score-recalibration.log
echo "Step 2 - Base Quality Score Recalibration - Finished"
echo "Step 2 outputs stored in: ${WKD}/base-quality-score-recalibration/"

# Step 3 - Collect alignment quality metrics
echo "Step 3 - Alignment quality - Started"
collect-alignment-metrics.sh \
    ${WKD}/base-quality-score-recalibration/${DATASET}_bqsr.bam \
    ${DATASET}_quality \
    ${FASTA} \
    ${SAMPLE_SIZE} \
    ${SCRIPTS} \
    2>&1 > ${WKD}/run-logs/collect-alignment-metrics.log
echo "Step 3 - Alignment quality - Finished"
echo "Step 3 outputs stored in: ${WKD}/alignment-quality-metrics/"

# Step 4 - Variant calling
echo "Step 4 - Variant calling - Started"
variant-calling.sh \
    ${WKD}/base-quality-score-recalibration/${DATASET}_bqsr.bam \
    ${DATASET} \
    ${FASTA} \
    ${INTERVALS} \
    2>&1 > ${WKD}/run-logs/variant-calling.log
echo "Step 4 - Variant calling - Finished"
echo "Step 4 outputs stored in: ${WKD}/variant-calling/"


# Step 5 - Variant filtering and quality metrics
echo "Step 5 - Variant filtering and quality metrics - Started"
variant-filtering.sh \
    ${WKD}/variant-calling/${DATASET}.vcf.gz \
    ${DATASET} \
    ${FASTA} \
    ${DBSNP} \
    ${SCRIPTS} \
    ${SAMPLE_SIZE} \
    2>&1 > ${WKD}/run-logs/variant-filtering.log
echo "Step 5 - Variant filtering and quality metrics - Finished"
echo "Step 5 outputs stored in: ${WKD}/variant-filtering/"

echo "Workflow finished successfully!

See HTML reports for
- BAM quality metrics at: ${WKD}/alignment-quality-metrics/
- variant quality metrics at: ${WKD}/variant-filtering/"
