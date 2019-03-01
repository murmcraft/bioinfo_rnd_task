#!/bin/env Rscript --no-save
library(rmarkdown)
library(kableExtra)

args <- commandArgs(TRUE)
filename <- args[1]

outfile <- paste(filename, Sys.Date(), "BAM_QC.html", sep=".")
render(input = "bam_QC_document.Rmd", output_file = outfile)