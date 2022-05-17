#!/usr/bin/R

# output csv with mean of each file (csv)  

args = commandArgs(trailingOnly=TRUE)
csv_path <- args[1]

name <- gsub(".csv","", csv_path)
df <- read.csv(csv_path)

p_mean = round(mean(df$P), digits = 3)
p_noovlp = round(mean(df$P_noovlp), digits = 3)
range_p = round(range(df$P), digits = 3)
ovlp_len = median(df$original_refgene) - median(df$len_noovlp)

row <- paste(name, p_mean, p_noovlp, range_p[1], range_p[2], range_p[2]-range_p[1], ovlp_len)

#Writing out the CSV (append)
out.file<- file("./all_pdist.csv", "a")
cat(row, file=out.file, append=TRUE, sep = "\n")
close(out.file)
