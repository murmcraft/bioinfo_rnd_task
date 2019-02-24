#!/bin/bash

# Mark duplicates
# Mark duplicate reads (PCR, optical) to be ignored in variant calling

# User input variables
INBAM=$1 # Input .bam filename
OUT_PREFIX=$2 # Output file prefix

# Run Picard MarkDuplicates to add a flag to
# primary alingments which are duplicates
gatk MarkDuplicates \
    -I ${INBAM} \
    -O ${OUT_PREFIX}.bam \
    --METRICS_FILE ${OUT_PREFIX}.txt
