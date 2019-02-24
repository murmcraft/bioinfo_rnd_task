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

# Picard's CollectWgsMetrics
gatk CollectWgsMetrics \
	-I ${INBAM} \
	-O ${OUT_PREFIX}_wgs.metrics \
	-R ${FASTA} \
	--INCLUDE_BQ_HISTOGRAM true \
	--READ_LENGTH ${READ_LENGTH} \
	--USE_FAST_ALGORITHM true

# Picard's CollectQualityScoreDistribution
gatk QualityScoreDistribution \
	-I ${INBAM} \
	-O ${OUT_PREFIX}_qscore.metrics \
	--CHART ${OUT_PREFIX}_qscore.pdf \
	--ALIGNED_READS_ONLY true \
	--PF_READS_ONLY true