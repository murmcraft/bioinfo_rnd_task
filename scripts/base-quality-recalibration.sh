#!/bin/bash

# Base quality recalibration

# User input variables
INBAM=$1 # Input .bam filename
OUT_PREFIX=$2 # Output file prefix
FASTA=$3 # Reference genome .fasta file
SITES=$4 # Known variant sites VCF
INDELS=$5 # Known indels VCF

# Run BaseRecalibrator to build a covariation model
gatk BaseRecalibrator \
    -I ${INBAM} \
    -O ${OUT_PREFIX}.table \
    -R ${FASTA} \
    --use-original-qualities \
    --known-sites ${SITES} \
    --known-sites ${INDELS}

# Apply the model to adjust base calling quality scores
gatk ApplyBQSR \
	-I ${INBAM} \
    -O ${OUT_PREFIX}.bam \
    -R ${FASTA} \
    --use-original-qualities \
    --bqsr-recal-file ${OUT_PREFIX}.table \
    --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30
