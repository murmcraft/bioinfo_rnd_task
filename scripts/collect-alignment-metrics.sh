#!/bin/bash

# Read and alignment quality summary metrics

# User input variables
INBAM=$1 # Input .bam filename
OUT_PREFIX=$2 # Output file prefix
FASTA=$3 # Reference genome .fasta file
READ_LENGTH=$4 # Average read length (bp)
QUALITIES=$5 # false if BaseRecalibration is applied and should use OQ tag, otherwise true

# Collect samtools stats table
samtools stats \
    -r ${FASTA} \
    --remove-dups \
    ${INBAM} \
    > ${OUT_PREFIX}_samtools.metrics

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
    grep "^SN" | cut -f1-3 > ${OUT_PREFIX}_overall.stats
cat ${OUT_PREFIX}_samtools.metrics | \
    grep "^GCC\|^FBC\|^LBC" | cut -f1-6 > ${OUT_PREFIX}_acgt.stats
cat ${OUT_PREFIX}_samtools.metrics | \
    grep "^ID" | cut -f1-4 > ${OUT_PREFIX}_indeldistribution.stats
cat ${OUT_PREFIX}_samtools.metrics | \
    grep "^COV" | cut -f1,3-4 > ${OUT_PREFIX}_coverage.stats