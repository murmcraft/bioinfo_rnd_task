#!/bin/env Rscript --no-save
library(rmarkdown)
library(kableExtra)

args <- commandArgs(TRUE)
output_dir <- args[1]
filename <- args[2]
scripts <- args[3]
samplesize <- args[4]

outfile <- paste0(output_dir, "/", filename, ".", 
	Sys.Date(), ".variant_QC.html")
render(input = paste0(scripts, "/variant-QC-report.Rmd"), 
       output_file = outfile,
       params = list(scripts = scripts, 
       	input_dir = output_dir, 
       	dataset = filename,
       	sample = samplesize))
