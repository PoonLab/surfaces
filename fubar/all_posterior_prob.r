#!/usr/bin Rscript --vanilla

# Load the package required to read JSON files.
library("rjson")
library("ggplot2")

# input arguments 
args = commandArgs(trailingOnly=TRUE)
json_path <- args[1]
outpath <- args[2]

# Give the input file name to the function.
result <- fromJSON(file = json_path)

#grid as df
grid.matix <- matrix(unlist(result$grid),ncol=3,byrow=TRUE)
grid.df <- as.data.frame(grid.matix)
colnames(grid.df) <- c("dN","dS","posterior")

#as factor 
grid.factor <- grid.df
grid.factor$dN<- as.factor(round(grid.factor$dN,digits = 2))
grid.factor$dS <- as.factor(round(grid.factor$dS,digits = 2))

grid.factor$log <- log10(grid.factor$posterior)

# writing out the data frame (appending each to the bottom)
write.table(grid.factor, outpath, sep = ",", append = T, col.names = FALSE, row.names = FALSE, quote = FALSE)
