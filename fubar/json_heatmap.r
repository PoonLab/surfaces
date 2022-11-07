#!/usr/bin/R
# Load the package required to read JSON files.
library("rjson")
library("ggplot2")

### INPUT 
midpoint =-5.14635
limits=c(-9.5864, -0.7063)

# log default : midpoint =-5.14635, limits=c(-9.5864, -0.7063))
# log 10 : midpoint =-2.235, limits=c(-4.1633, -0.3067))

# 2. name_cleanup (gsub)
name <- gsub(".clean.FUBAR.json","",json_path)

# 3. file in and out
args = commandArgs(trailingOnly=TRUE)
json_path <- args[1]
outpath <- args[2]

# Give the input file name to the function.
result <- fromJSON(file = json_path)

#grid as df
grid.matix <- matrix(unlist(result$grid),ncol=3,byrow=TRUE)
grid.df <- as.data.frame(grid.matix)
### colnames(grid.df) <- c("dN","dS","posterior") # according to the documentation
colnames(grid.df) <- c("dS","dN","posterior") # confirmed this is the right lables with Art 

#as factor 
grid.factor <- grid.df
grid.factor$dN<- as.factor(round(grid.factor$dN,digits = 3))
grid.factor$dS <- as.factor(round(grid.factor$dS,digits = 3))
grid.factor$log <- log(grid.factor$posterior) # plot log transformed posterior probability base e
grid.factor$log10 <- log(grid.factor$posterior,10) # plot log transformed posterior probability base 10
grid.df$name <- name

# Heatmap 
png(file=outpath,width=800, height=700)
ggplot(grid.factor, aes(grid.factor$dS, grid.factor$dN, fill= grid.factor$log)) +
geom_tile() + ggtitle(name) +
  xlab("Synonymous rate (dS)") + 
ylab("Nonsynonymous rate (dN)") + 
labs(fill = "Posterior Mean") +
theme(axis.text.x = element_text(angle = 90), 
      axis.text = element_text(size=13),
      axis.title = element_text(size=16,face="bold"),
      plot.title = element_text(lineheight=2, 
                                face="bold", size = 20, hjust = 0.5)) + 
scale_fill_gradient2('posterior', low = "blue", mid = "white", high = "red", 
                     midpoint = mid_range, limits= prob_range + 
                     coord_fixed() +
                     guides(fill = guide_colourbar(barwidth = 2,
                                barheight = 40,
                                title = "Posterior Mean",
                                ticks = FALSE))
dev.off()




"""
# Writing it out to a file all_posterior_prob.csv
#output.file <- file("./all_posterior-prob.csv", "a")
#grid.df$log <- log(grid.df$posterior)
#colnames(grid.df) <- c("dN","dS","posterior","log")
#write.table(grid.df, output.file, sep = ",", append = T)
#close(output.file)

# max & min 
#max_min <- paste(name,max(grid.df$posterior),min(grid.df$posterior))
#max.output.file <- file("./max.min.pearson_cor.csv", "a")
#cat(max_min, file=max.output.file, append=TRUE, sep = "\n")
#close(max.output.file)
"""
                                          
                              
### INFO
# Read jason file and save jpg of heatmap of posterior probability with set range
# old_ncbi: midpoint = -5.175, limits=c(-9.59, -0.76))
# lin(april): midpoint = -2.12, limits=c(-3.86, -0.37))
# ncbi(april): midpoint =-2.06, limits=c(-3.80, -0.31))
                     
                     
