#!/usr/bin Rscript --vanilla
# Load the package required to read JSON files.

library(seqinr)
require(seqinr)

args = commandArgs(trailingOnly=TRUE)
msa.path <- args[1]

#clean the file name to only be virus_accn_gene_accn (ncbi)
name <- gsub("cutter_cds_", "", msa.path)
name <- gsub(".clean", "", msa.path)

mx <- do.call(rbind, aln) # as a matrix
vec <- as.vector(a) # as a vector 

# Total num of gaps in all seq
num_gap <- length(which(vec=="-"))
  
# Total number of nuc in all sequences 
total_length <- length(vec)

#gap percentage 
gap_per <- num_gap/total_length

# mean gap percentage of all sequences of a cds 
result <- paste(name, gap_per)

# Writing it out to a file 
output.file <- file("./gap_per.csv", "a")
cat(result, file=output.file, append=TRUE, sep = "\n")
close(output.file)
print(result)
