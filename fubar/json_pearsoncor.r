###################################### ART REVISED VERSION ########
#!/usr/bin Rscript --vanilla
# Load the package required to read JSON files.

library("rjson")
setwd(dir="/home/sareh/re-do/FUBAR/fubar100/json/")

args = commandArgs(trailingOnly=TRUE)
json_path1 <- args[1]
json_path2 <- args[2]

import.json <- function(path) {
  result <- fromJSON(file = path)
  matix <- matrix(unlist(result$grid), ncol=3, byrow=TRUE)
  as.data.frame(matix)
}

grid.df1 <- import.json(json_path1)
grid.df2 <- import.json(json_path2)

name1 <- gsub("clean.FUBAR.json", "", json_path1)
name2 <- gsub("clean.FUBAR.json", "", json_path2)

# quality control
stopifnot(dim(grid.df1) == dim(grid.df2))
stopifnot(all(grid.df1$V1 == grid.df2$V1))
stopifnot(all(grid.df1$V2 == grid.df2$V2))

# Pearson Correlation 
res <- cor.test(grid.df1$V3, grid.df2$V3, method = "pearson")

result = paste(name1, name2, 1-res$estimate)

# Writing it out to a file 
output.file <- file("./pearson_correlation.csv", "a")
cat(result, file=output.file, append=TRUE, sep = "\n")
close(output.file)
print(result)
