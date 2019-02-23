#!/bin/bash

# Read and alignment quality summary metrics

# User input variables
INBAM=$1 # Input .bam filename
OUT_PREFIX=$2 # Output file prefix
FASTA=$3 # Reference genome .fasta file
READ_LENGTH=$4 # Average read length (bp)

# First, run Picard's CollectAlignmentSummaryMetrics
gatk CollectAlignmentSummaryMetrics \
	-I ${INBAM} \
	-O ${OUT_PREFIX}_calnm.txt \
	-R ${FASTA}

# Second, run Picard's CollectWgsMetrics
gatk CollectWgsMetrics \
	-I ${INBAM} \
	-O ${OUT_PREFIX}_cwgsm.txt \
	-R ${FASTA} \
	--INCLUDE_BQ_HISTOGRAM true \
	--READ_LENGTH ${READ_LENGTH} \
	--USE_FAST_ALGORITHM true
