#!/bin/bash

# Preparatory step - File validation
# The input BAM file must be validated before running the pipeline

# User input variables
INBAM=$1 # Path to input .bam file
OUT_PREFIX=$2 # Output file prefix

# Create an output directory
DIR=validate-sam-file 
rm -rf $DIR 
mkdir -p $DIR 
cd $DIR

# Run Picard ValidateSamFile to diagnose potential issues
# First see summary
gatk ValidateSamFile \
    -I ${INBAM} \
    -MODE SUMMARY | \
grep "ERROR:" \
> ${OUT_PREFIX}.summary

# In case there were errors, store the details into a txt file
if [[ -s ${OUT_PREFIX}.summary ]]; then
    gatk ValidateSamFile \
        -I ${INBAM} \
        -IGNORE_WARNINGS true \
        -MODE VERBOSE \
    > ${OUT_PREFIX}.verbose
else 
    rm ${OUT_PREFIX}.summary
fi

# Fix the identified errors, GATK documentation points some directions:
# https://software.broadinstitute.org/gatk/documentation/article.php?id=7571