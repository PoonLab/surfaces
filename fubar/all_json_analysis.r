#!/usr/bin Rscript --vanilla

# input : json 
# output 1. all_posterior_prob.csv 
# output 2. is a file of MLE.df [outpath]
# output 3. is appended to a csv [all_mle.csv] 

# Load the package required to read JSON files.
library("ggplot2")

# opening json with jsonlite
library(jsonlite)
result <- read_json(json_path)

# INPUT 
name <- gsub(".FUBAR.json","", json_path) #file name

# log transformed ranges
args = commandArgs(trailingOnly=TRUE)
json_path <- args[1]
outpath <- args[2]


### MLE ###
MLE.matix <- matrix(unlist(result$MLE$content),ncol=8,byrow=TRUE)
MLE.df <- as.data.frame(MLE.matix)
MLE.df <- MLE.df[,c(1:6)]
MLE.df$beta.alpha <- MLE.df$beta/MLE.df$alpha # new
colnames(MLE.df) <- c("alpha","beta","beta-alpha","Prob[alpha>beta]","Prob[alpha<beta]","BayesFactor[alpha<beta]","beta/alpha")          

#info
codon <- nrow(MLE.df)
mean.dS <- mean(MLE.df$alpha)
mean.dN <- mean(MLE.df$beta)
pos.sites <- nrow(MLE.df[MLE.df$`Prob[alpha<beta]`>0.8,])
neg.sites <- nrow(MLE.df[MLE.df$`Prob[alpha>beta]`>0.8,])
mean.beta-alpha <- mean(MLE.df$beta-alpha) #new
mean.beta/alpha <- mean(MLE.df$beta/alpha) #new 

row <- paste(name, codon, mean.dN, mean.dS, pos.sites, neg.sites, mean.dN/mean.dS, pos.sites/codon,neg.sites/codon,mean.beta-alpha,mean.beta/alpha)                      
print(row)

### GRID ###
#grid as df
grid.matix <- matrix(unlist(result$grid),ncol=3,byrow=TRUE)
grid.df <- as.data.frame(grid.matix)
colnames(grid.df) <- c("dN","dS","posterior")
grid.df$log <- log(grid.df$posterior)

#as factor 
grid.factor <- grid.df
grid.factor$dN<- as.factor(round(grid.factor$dN,digits = 2))
grid.factor$dS <- as.factor(round(grid.factor$dS,digits = 2))
grid.factor$log <- log10(grid.factor$posterior)



### OUTPUT FILES ###
# 1.all_posterior_prob.csv: writing out the data frame (appending each to the bottom)
write.table(grid.factor, outpath, sep = ",", append = T)

# Writing it out MLE
output.file <- file(outpath, "w")
write.table(MLE.df, output.file, sep = ",")
close(output.file)

#Writing out the CSV (append)
out.file<- file("./all_mle.csv", "a")
cat(row, file=out.file, append=TRUE, sep = "\n")
close(out.file)

# Heatmap auto range 
png(file=outpath,width=800, height=700)
ggplot(grid.factor, aes(grid.factor$dS, grid.factor$dN, fill= grid.factor$log)) +
geom_tile() + theme(axis.text.x = element_text(angle = 90),
axis.text=element_text(size=13),axis.title=element_text(size=16,face="bold")) +
xlab("Synonymous rate (dS)") +
ylab("Nonsynonymous rate (dN)") +
labs(fill = "Posterior Mean") +
ggtitle(name) + theme(plot.title = element_text(lineheight=2, face="bold", size = 20, hjust = 0.5)) +
scale_fill_gradient2('posterior', low = "blue", mid = "white", high = "red") +
coord_fixed() +
guides(fill = guide_colourbar(barwidth = 2,barheight = 40,title = "Posterior Mean",ticks = FALSE))
dev.off()

# max & min 
max_min <- paste(name, max(grid.df$posterior), min(grid.df$posterior), max(grid.df$log), min(grid.df$log))
max.output.file <- file("./max.min.pearson_cor.csv", "a")
cat(max_min, file=max.output.file, append=TRUE, sep = "\n")
close(max.output.file)


# Heatmap constant range 

mid_range = -2.2454638
prob_range = c(-4.1639958,-0.3269318)

png(file=outpath,width=800, height=700)
ggplot(grid.factor, aes(grid.factor$dS, grid.factor$dN, fill= grid.factor$log)) +
geom_tile() + theme(axis.text.x = element_text(angle = 90),
axis.text=element_text(size=13),axis.title=element_text(size=16,face="bold")) +
xlab("Synonymous rate (dS)") +
ylab("Nonsynonymous rate (dN)") +
labs(fill = "Posterior Mean") +
ggtitle(name) + theme(plot.title = element_text(lineheight=2, face="bold", size = 20, hjust = 0.5)) +
scale_fill_gradient2('posterior', low = "blue", mid = "white", high = "red", midpoint = mid_range, limits = prob_range) +
coord_fixed() +
guides(fill = guide_colourbar(barwidth = 2,barheight = 40,title = "Posterior Mean",ticks = FALSE))
dev.off()



