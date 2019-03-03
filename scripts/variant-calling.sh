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
CPU_COUNT=2 # Hard-coded for 8 CPU laptop
PID_COUNTER=0 
cat ${INTERVAL_LIST} | \
while read -r INTERVAL; do
    PID_COUNTER=$((${PID_COUNTER} + 1))
    # Run the variant calling in background
    (gatk HaplotypeCaller \
        -I ${INBAM} \
        -O ${OUT_PREFIX}_${INTERVAL}.vcf.gz \
        -R ${FASTA} \
        -L ${INTERVAL};
        # In case the run fails, collect the interval to a fail list
        if [[ $? -ne 0 ]]; then 
            echo "${INTERVAL}" >> interval.fail
        fi
    ) &
    # Run only 2 intervals simultaneously in background
    if [[ ${PID_COUNTER} -gt ${CPU_COUNT} ]]; then
        wait
        PID_COUNTER=0
    fi
done
# Wait for all the interval processes to finish
wait

# See if any of the interval runs failed, if yes, stop the run
if [[ -s interval.fail ]]; then
    exit 1
fi

# Concatenate all the intervals into a single VCF
ls *.vcf.gz | sort -V > interval_vcfs.list
bcftools concat -f interval_vcfs.list \
    -Oz -o ${OUT_PREFIX}.vcf.gz
# Create index
bcftools index -t ${OUT_PREFIX}.vcf.gz

# Clean up
rm ${OUT_PREFIX}_*:*.vcf.gz*
rm interval_vcfs.list