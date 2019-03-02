#!/bin/bash

# Variant calling 
# Variant calling with GATK HaplotypeCaller per genomic interval
# and finally merging all individual files together

# User input variables
INBAM=$1 # Path to input .bam file
OUT_PREFIX=$2 # Output file prefix
FASTA=$3 # Reference genome .fasta file
INTERVAL_LIST=$4 # Interval file containing whole genome good intervals

# Create an output directory
DIR=variant-calling 
rm -rf $DIR 
mkdir -p $DIR
cd $DIR

# Variant calling per interval
CPU_COUNT=4
PID_COUNTER=0
cat ${INTERVAL_LIST} | \
while read -r INTERVAL; do
    PID_COUNTER=$((${PID_COUNTER} + 1))
    (gatk HaplotypeCaller \
        -I ${INBAM} \
        -O ${OUT_PREFIX}_${INTERVAL}.vcf.gz \
        -R ${FASTA} \
        -L ${INTERVAL};
        if [[ $? -ne 0 ]]; then 
            echo "${INTERVAL}" >> interval.fail
        fi
    ) &
    if [[ ${PID_COUNTER} -gt ${CPU_COUNT} ]]; then
        wait
        PID_COUNTER=0
    fi
done
wait

if [[ -s interval.fail ]]; then
    exit 1
fi