#!/usr/bin Rscript --vanilla
# Load the package required to read JSON files.

library(seqinr)
require(seqinr)

#setwd(dir="/home/sareh/new_cds/hyphyclean_cds")

args = commandArgs(trailingOnly=TRUE)
msa.path <- args[1]

aln <- read.fasta(msa.path, seqtype="AA")
mx <- do.call(rbind, aln)

ent <- function(x) {
  p <- table(as.character(x))
  p <- p[names(p) != '-']  # ignore gaps
  p <- p/sum(p)
  sum(-p*log(p))
}

ent.prof <- apply(mx, 2, ent)
mean <- mean(ent.prof)

name <- gsub("cutter_cds_", "", msa.path)
name <- gsub(".clean", "", msa.path)

result <- paste(name, mean)

# Writing it out to a file 
output.file <- file("./aa_entropy.csv", "a")
cat(result, file=output.file, append=TRUE, sep = "\n")
close(output.file)
print(result)
