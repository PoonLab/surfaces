#!/usr/bin/Rscript
# Read jason file and save jpg of heatmap of posterior probability 

#this didn't work: Error in file(con, "r") : cannot open the connection
args = commandArgs(trailingOnly=TRUE)

json_path <- args[1]
outpath <- args[2]

#title
name <- gsub("/home/sareh/re-do/FUBAR/fubar100/json/","",json_path)
name <- gsub("clean.FUBAR.json","",name)

# Load the package required to read JSON files.
library("rjson")
library("ggplot2")

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

# Heatmap 
png(file=outpath,width=600, height=600)
ggplot(grid.factor, aes(grid.factor$dS, grid.factor$dN, fill= grid.factor$posterior)) +
geom_tile() + theme(axis.text.x = element_text(angle = 90)) +
  xlab("Synonymous rate") + ylab("Nonsynonymous rate") + labs(fill = "Posterior Mean") +
ggtitle(name) + theme(plot.title = element_text(lineheight=2, face="bold", size = 10, hjust = 0.5)
                      + scale_fill_gradient(low = "blue", high = "red")
dev.off()

