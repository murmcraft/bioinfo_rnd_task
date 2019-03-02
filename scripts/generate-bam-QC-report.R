#!/bin/env Rscript --no-save
library(rmarkdown)
library(kableExtra)

args <- commandArgs(TRUE)
output_dir <- args[1]
filename <- args[2]
scripts <- args[3]

outfile <- paste0(output_dir, "/", filename, ".", 
	Sys.Date(), ".BAM_QC.html")
render(input = paste0(scripts, "/bam-QC-report.Rmd"), 
       output_file = outfile,
       params = list(scripts = scripts, 
       	input_dir = output_dir))
