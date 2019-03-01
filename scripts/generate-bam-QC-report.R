#!/bin/env Rscript --no-save
library(rmarkdown)
library(kableExtra)

args <- commandArgs(TRUE)
filename <- args[1]

outfile <- paste(filename, Sys.Date(), "BAM_QC.html", sep=".")
render(input = "bam-QC-report.Rmd", output_file = outfile)