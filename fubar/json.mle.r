#!/usr/bin/R

# output 1. is a file of MLE.df [outpath]
# output 2. is appended to a csv [all_mle.csv] 

args = commandArgs(trailingOnly=TRUE)

json_path <- args[1]
outpath <- args[2]

#title
#cutter_cds_NC_039223_YP_009513196.1.clean.FUBAR.json
name <- gsub(".noovlp.clean.FUBAR.json","")

# Load the package required to read JSON files.
library("rjson")
library("ggplot2")

# Give the input file name to the function.
result <- fromJSON(file = json_path)

#MLE.matix
MLE.matix <- matrix(unlist(result$MLE$content),ncol=8,byrow=TRUE)
#MLE.data frame
MLE.df <- as.data.frame(MLE.matix)
MLE.df <- MLE.df[,c(1:6)]
colnames(MLE.df) <- c("alpha","beta","beta-alpha","Prob[alpha>beta]","Prob[alpha<beta]","BayesFactor[alpha<beta]")

codon <- nrow(MLE.df)
mean.dS <- mean(MLE.df$alpha)
mean.dN <- mean(MLE.df$beta)
pos.sites <- nrow(MLE.df[MLE.df$`Prob[alpha<beta]`>0.8,])
neg.sites <- nrow(MLE.df[MLE.df$`Prob[alpha>beta]`>0.8,])

row <- paste(name, codon, mean.dN, mean.dS, pos.sites, neg.sites, mean.dN/mean.dS, pos.sites/codon,neg.sites/codon)
print(row) 

# Writing it out MLE
output.file <- file(outpath, "w")
write.table(MLE.df, output.file, sep = ",")
close(output.file)

#Writing out the CSV (append)
out.file<- file("./all_mle.csv", "a")
cat(row, file=out.file, append=TRUE, sep = "\n")
close(out.file)
