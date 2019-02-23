# Tiny pipeline for simple variant calling and quality control

This pipeline implements the following:
- input as BAM file
- BAM quality control
- variant calling on split BAM files
- merging the variant results together
- variant calling summary report

The workflow is tested with a small WGS dataset - NA12878 aligned to GRCh37 reference genome.  
The raw data was obtained from [GATK test data sets](https://console.cloud.google.com/storage/browser/gatk-test-data).

## Build and run the Docker image

The Dockerfile is based on [Broad Institute's GATK docker](https://hub.docker.com/r/broadinstitute/gatk/) with added tools and functionalities.

It includes GATK4 and additional tools required for the pipeline:
- GATK v4.1.0.0
- SAMtools v1.9
- BCFtools v1.9
- R v3.2.5

First, clone the repository and enter the directory.

Build the image:
```
docker build -t bioinfo-rnd-task
```

Then, run it while mounting your input data directory (here `testdata`) to `/home`:
```
docker run -ti \
    -v /path/to/testdata/:/home \
    bioinfo-rnd-task
```

## Input file requirements

The input BAM file should be:
- paired-end WGS data,
- sorted,
- indexed, and
- in correct format according to GATK requirements.

To confirm that your data is indeed in suitable format,  
**run the validation and fix errors**, if necessary:
```
validate-sam-file.sh \
    NA12878.bam \
    validate
```
where `NA12878.bam` is input BAM filename and `validate` is output filename prefix.

In case errors where found, the script procudes `.summary`and `.verbose` files, which indicate the errors and point their sources.  
[GATK documentation](https://software.broadinstitute.org/gatk/documentation/article.php?id=7571) provides tips how to fix the BAM file into a compatible format. 


## Run the full workflow

To run the entire workflow, use the script with wanted inputs:

*********UPDATE HERE

The full workflow will run each step described below.  
Wanted steps can also be run separately according to the command examples.

## Read and alignment quality



First, collect alignment quality metrics:
```
collect-alignment-metrics.sh \
    NA12878.bam \
    alignment-metrics \
    Homo_sapiens.GRCh37.dna.primary_assembly.fa \
    250
```
where `NA12878.bam` is input BAM filename, `alignment-metrics` is output filename prefix, `Homo_sapiens.GRCh37.dna.primary_assembly.fa` is the reference genome fasta file, and `250` is the average read length.

