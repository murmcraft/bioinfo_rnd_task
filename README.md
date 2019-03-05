# Tiny pipeline for simple variant calling and quality control

This pipeline implements the following:
- input as BAM file
- BAM quality control
- variant calling on split BAM files
- merging the variant results together
- variant calling summary report

The workflow is tested with a small WGS dataset - NA12878 aligned to GRCh37 reference genome.  
The raw data was obtained from [GATK test data sets](https://console.cloud.google.com/storage/browser/gatk-test-data) and aligned to the GRCh37 reference genome with `bwa mem` according to [GATK workflow](https://github.com/gatk-workflows/gatk4-data-processing/blob/master/processing-for-variant-discovery-gatk4.wdl).

## Build and run the Docker image

The Dockerfile is based on [Broad Institute's GATK docker](https://hub.docker.com/r/broadinstitute/gatk/) with added tools and functionalities.

It e.g. includes GATK4 and additional tools required for the pipeline:
- GATK v4.1.0.0
- SAMtools v1.9
- BCFtools v1.9
- R v3.2.5

First, clone the repository and enter the directory.

Build the image:
```
docker build -t bioinfo-rnd-task .
```

Then, run it while mounting your input data directory (here `testdata`) to `/home`:
```
docker run -ti \
    -v /path/to/testdata/:/home \
    bioinfo-rnd-task
```

First, setup your working directory and ensure your input files are there:
```
WKD=/home/
ls ${WKD}
```

## Input file requirements and preparations

#### BAM file
The input BAM file should *at least* be:
- paired-end WGS data,
- trimmed or clipped for adapters,
- sorted,
- indexed, and
- contain read group(s) according to GATK requirements.

To confirm that your data indeed is in suitable format, 
**run the validation and fix errors**, if necessary:
```
validate-sam-file.sh \
    ${WKD}/NA12878.bam \
    NA12878_validate
```
where `/home/NA12878.bam` is path to input BAM file and `NA12878_validate` is output filename prefix.

In case errors were found, the script procudes 
`.summary`and `.verbose` files, which indicate the errors and point their sources.  
[GATK documentation](https://software.broadinstitute.org/gatk/documentation/article.php?id=7571) provides tips how to fix the BAM file into a compatible format. 

#### Reference files
The required files are:
- reference genome fasta with dictionary and index,
- known variant sites,
- known indels,
- truth datasets, and 
- WGS interval list

Depending on to which reference genome your reads are aligned, download corresponding reference files.  
Here, the test data is aligned to GRCh37 reference genome and the corresponding files can be downloaded from indicated sources:
```
mkdir -p reference-data
cd reference-data
wget ftp://ftp.ensembl.org/pub/grch37/release-95/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz
wget --user='gsapubftp-anonymous' --password='' \
    ftp://ftp.broadinstitute.org/bundle/b37/dbsnp_138.b37.excluding_sites_after_129.vcf.gz*
wget --user='gsapubftp-anonymous' --password='' \
    ftp://ftp.broadinstitute.org/bundle/b37/Mills_and_1000G_gold_standard.indels.b37.vcf.gz*
wget --user='gsapubftp-anonymous' --password='' \
    ftp://ftp.broadinstitute.org/bundle/b37/1000G_omni2.5.b37.vcf.gz*
wget --user='gsapubftp-anonymous' --password='' \
    ftp://ftp.broadinstitute.org/bundle/b37/1000G_phase3_v4_20130502.sites.vcf.gz*
```
Unzip, BGZF compress and index the files:
```
gunzip Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz
gunzip dbsnp_138.b37.excluding_sites_after_129.vcf.gz
gunzip Mills_and_1000G_gold_standard.indels.b37.vcf.gz
gunzip 1000G_omni2.5.b37.vcf.gz
gunzip 1000G_phase3_v4_20130502.sites.vcf.gz
bcftools view dbsnp_138.b37.excluding_sites_after_129.vcf \
    -Oz -o dbsnp_138.b37.excluding_sites_after_129.vcf.gz
bcftools index -t dbsnp_138.b37.excluding_sites_after_129.vcf.gz
bcftools view Mills_and_1000G_gold_standard.indels.b37.vcf \
    -Oz -o Mills_and_1000G_gold_standard.indels.b37.vcf.gz
bcftools index -t Mills_and_1000G_gold_standard.indels.b37.vcf.gz
bcftools view 1000G_omni2.5.b37.vcf -Oz -o 1000G_omni2.5.b37.vcf.gz
bcftools index -t 1000G_omni2.5.b37.vcf.gz
bcftools view 1000G_phase3_v4_20130502.sites.vcf \
    -Oz -o 1000G_phase3_v4_20130502.sites.vcf.gz
bcftools index -t 1000G_phase3_v4_20130502.sites.vcf.gz
```
To generate the fasta file dictionary, use:
```
gatk CreateSequenceDictionary \
    -R Homo_sapiens.GRCh37.dna.primary_assembly.fa
    -O Homo_sapiens.GRCh37.dna.primary_assembly.dict
```
To generate the fasta file index, use:
```
samtools faidx Homo_sapiens.GRCh37.dna.primary_assembly.fa
```
Download the WGS interval list manually from [GATK Google Cloud bucket](https://console.cloud.google.com/storage/browser/gatk-test-data/intervals) and place it to the reference-data directory
- `b37_wgs_consolidated_calling_intervals.list`

## Run the complete workflow

To run the complete workflow, use the script with wanted inputs:

*********UPDATE HERE

The full workflow will run each step described below. Wanted steps can also be run separately according to the command examples.


## Step-by-step descriptions

### Mark duplicates

It is recommended practice to mark duplicate reads, which are then ignored in variant calling. Duplicates are determined as those reads whose 5' positions are identical. Often the origin of duplicate reads is the PCR during library preparation, but also sequencing may cause optical duplicates.

Run the script to flag duplicate reads and generate a summary metrics file:
```
mark-duplicates.sh \
    ${WKD}/NA12878.bam \
    NA12878_markdups
```
where `${WKD}/NA12878.bam` is the path to input BAM file and `NA12878_markdups` is a prefix for the output `.bam` and `.txt` files. 


### Base quality score recalibration

Base quality score recalibration is meant to detect systematic errors in the base calling quality scores made by the sequencer. The quality scores are adjusted with the help of known variants. 

```
base-quality-score-recalibration.sh \
    ${WKD}/mark-duplicates/NA12878_markdups.bam \
    NA12878_bqsr \
    /home/reference-data/Homo_sapiens.GRCh37.dna.primary_assembly.fa \
    /home/reference-data/dbsnp_138.b37.excluding_sites_after_129.vcf.gz \
    /home/reference-data/Mills_and_1000G_gold_standard.indels.b37.vcf.gz
```
where `${WKD}/NA12878.bam` is the path to input BAM file, `NA12878_bqsr` is a output filename prefix, `Homo_sapiens.GRCh37.dna.primary_assembly.fa` is the reference genome fasta file, `dbsnp_138.b37.excluding_sites_after_129.vcf.gz` is the dbSNP known variants and `Mills_and_1000G_gold_standard.indels.b37.vcf.gz` contain the known indels. 


### Alignment quality

Collect various alignment quality metrics:
```
collect-alignment-metrics.sh \
    ${WKD}/base-quality-score-recalibration/NA12878_bqsr.bam \
    NA12878_quality \
    /home/reference-data/Homo_sapiens.GRCh37.dna.primary_assembly.fa \
    0.99 \
    /scripts
```
where `NA12878_bqsr.bam` is input BAM filename, `NA12878_quality` is output filename prefix, `Homo_sapiens.GRCh37.dna.primary_assembly.fa` is the reference genome fasta file, `0.99` is the subsampling proportion (here 99 %, because the test dataset is extremely small, but for a full WGS dataset one should choose a small value e.g. 0.25 depending on the size of the input data), and `/scripts` is the path to scripts directory.

The script generates an HTML report `NA12878_bqsr.<date>.BAM_QC_report.html` containing BAM quality metrics (currently only summary metrics table, ACGT content per cycle plot, coverage, mapping quality and indel lengths histograms are implemented).


### Variant calling

GATK HaplotypeCaller simultaneously calls SNPs and indels and does local *de novo* assembly at the active region increasing the accuracy of the calls. To speed up the variant calling, it is good to parallelize the process per genomic interval. These are defined by GATK as a handful of non-overlapping regions per each chromosome excluding non-interesting or difficult regions such as centromeres (in total 103 intervals of autosomal and X chromosomes).  

Here, due to the tiny example run on a laptop (8 logical cores), the parallelization is hard-coded to 2 simultaneous processes, where GATK by default uses 4 threads for each process. For a real data run for instance on a cluster, all 103 jobs could be simultaneously submitted to a workload manager. 

Call variants per genomic intervals:
```
variant-calling.sh \
    ${WKD}/base-quality-score-recalibration/NA12878_bqsr.bam \
    NA12878 \
    /home/reference-data/Homo_sapiens.GRCh37.dna.primary_assembly.fa \
    /home/reference-data/intervals_b37_wgs_consolidated_calling_intervals.list
```

