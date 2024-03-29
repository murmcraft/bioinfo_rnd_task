#!/bin/bash

# Read and alignment quality summary metrics

# User input variables
INBAM=$1 # Input .bam filename
OUT_PREFIX=$2 # Output file prefix
FASTA=$3 # Reference genome .fasta file
SAMPLE_SIZE=$4 # As fraction for downsampling a huge input file
SCRIPTS=$5 # Path to scripts

# Create an output directory
DIR=alignment-quality-metrics 
rm -rf $DIR 
mkdir -p $DIR 
cd $DIR

# Collect samtools stats table
samtools stats \
    -r ${FASTA} \
    --remove-dups \
    ${INBAM} \
    > ${OUT_PREFIX}_samtools.metrics

# Mapping quality distribution for a subsample
samtools view -s ${SAMPLE_SIZE} ${INBAM} | \
    cut -f5 \
    > ${OUT_PREFIX}_mapq.stats

# Picard's CollectAlignmentSummaryMetrics
gatk CollectAlignmentSummaryMetrics \
    -I ${INBAM} \
    -O ${OUT_PREFIX}_alignment.metrics \
    -R ${FASTA}

# Picard's CollectQualityScoreDistribution
gatk QualityScoreDistribution \
    -I ${INBAM} \
    -O ${OUT_PREFIX}_qscore.metrics \
    --CHART ${OUT_PREFIX}_qscore.pdf \
    --ALIGNED_READS_ONLY true \
    --PF_READS_ONLY true

# Subset the reports into plottable formats
cat ${OUT_PREFIX}_samtools.metrics | \
    grep "^GCC\|^FBC\|^LBC" | cut -f1-6 > ${OUT_PREFIX}_acgt.stats
cat ${OUT_PREFIX}_samtools.metrics | \
    grep "^ID" | cut -f1-4 > ${OUT_PREFIX}_indeldistribution.stats
cat ${OUT_PREFIX}_samtools.metrics | \
    grep "^COV" | cut -f1,3-4 > ${OUT_PREFIX}_coverage.stats
cat ${OUT_PREFIX}_qscore.metrics | \
    grep '^[0-9]' > ${OUT_PREFIX}_qscore.stats

# Generate a HTML report of the quality metrics
FILE_PREFIX=$(basename ${INBAM%.bam})
OUTPUT_DIR=$(pwd)
Rscript ${SCRIPTS}/generate-bam-QC-report.R \
    ${OUTPUT_DIR} \
    ${FILE_PREFIX} \
    ${SCRIPTS}

# Clean up
rm ${OUT_PREFIX}_qscore.pdf
rm ${OUT_PREFIX}*.stats