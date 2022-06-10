#!/usr/bin/R

### INPUT 
# json
# midpoint =-2.06, limits=c(-3.80, -0.31))
# limits
# name_cleanup (gsub)

### OUTPUT
# heatmap.png
# all_posterior-prob.csv
# max.min.pearson_cor.csv

### INFO
# Read jason file and save jpg of heatmap of posterior probability with set range
# old_ncbi: midpoint = -5.175, limits=c(-9.59, -0.76))
# lin(april): midpoint = -2.12, limits=c(-3.86, -0.37))
# ncbi(april): midpoint =-2.06, limits=c(-3.80, -0.31))

args = commandArgs(trailingOnly=TRUE)
json_path <- args[1]
outpath <- args[2]

# log transformed ranges
mid_range = -2.2454638
prob_range = c(-4.1639958,-0.3269318)

#title
#cutter_cds_NC_039223_YP_009513196.1.clean.FUBAR.json
name <- gsub(".clean.FUBAR.json","",gsub("/home/sareh/new_cds/FUBAR/json/cutter_cds_","",json_path))

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
grid.factor$log <- log(grid.factor$posterior) # plot log transformed posterior probability
grid.df$name <- name


# Heatmap 
png(file=outpath,width=800, height=700)
ggplot(grid.factor, aes(grid.factor$dS, grid.factor$dN, fill= grid.factor$log)) +
geom_tile() + theme(axis.text.x = element_text(angle = 90), axis.text=element_text(size=13),axis.title=element_text(size=16,face="bold")) +
  xlab("Synonymous rate (dS)") + ylab("Nonsynonymous rate (dN)") + labs(fill = "Posterior Mean") +
ggtitle(name) + theme(plot.title = element_text(lineheight=2, face="bold", size = 20, hjust = 0.5)) +  
  
scale_fill_gradient2('posterior', low = "blue", mid = "white", high = "red",midpoint = mid_range, limits= prob_range + coord_fixed() +
  
  guides(fill = guide_colourbar(barwidth = 2,
                                barheight = 40,
                                title = "Posterior Mean",
                                ticks = FALSE))

dev.off()


"""
# Writing it out to a file all_posterior_prob.csv
output.file <- file("./all_posterior-prob.csv", "a")
grid.df$log <- log(grid.df$posterior)
colnames(grid.df) <- c("dN","dS","posterior","log")
write.table(grid.df, output.file, sep = ",", append = T)
close(output.file)

# max & min 
max_min <- paste(name,max(grid.df$posterior),min(grid.df$posterior))
max.output.file <- file("./max.min.pearson_cor.csv", "a")
cat(max_min, file=max.output.file, append=TRUE, sep = "\n")
close(max.output.file)
"""
