#!/bin/bash

# Mark duplicates
# Mark duplicate reads (PCR, optical) to be ignored in variant calling

# User input variables
INBAM=$1 # Path to input .bam file
OUT_PREFIX=$2 # Output file prefix

# Create an output directory
DIR=mark-duplicates 
rm -rf $DIR 
mkdir -p $DIR
cd $DIR

# Run Picard MarkDuplicates to add a flag to
# primary alingments which are duplicates
gatk MarkDuplicates \
    -I ${INBAM} \
    -O ${OUT_PREFIX}.bam \
    --METRICS_FILE ${OUT_PREFIX}.txt
