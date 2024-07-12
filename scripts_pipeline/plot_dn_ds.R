# Plot dN (beta) and dS (alpha) from FUBAR analysis
# Analyse FUBAR results for multiple CDSs in a virus
library(plyr)
require(jsonlite)
library(dplyr)

##################################
# Import JSON files 
##################################
setwd("/home/laura/Projects/surfaces/scripts_pipeline/temp_zika/zika_results")
files <- Sys.glob("*.FUBAR.json")
names(files) <- files  # so master list has filenames as names
master <- lapply(files, function(f) read_json(f, simplifyVector=TRUE))
names(master)

# Returns MLE estimates for each dataset as list of data frames
# Data frame res[[1]] includes alpha, beta, alpha=beta, LRT, P-value, Total branch length, codon.pos and dnds
res <- lapply(master, function(d) {
  stats <- d$MLE
  site.stats <- as.data.frame(stats$content[[1]][,1:6])
  names(site.stats) <- stats$headers[,1]
  site.stats$codon.pos<-row.names(site.stats)
  site.stats$dnds<- site.stats$beta/site.stats$alpha
  return(site.stats)
})

###############################################
# distribution of dN and dS for all proteins
###############################################
par(mfrow=c(4,4)) # 6 rows, 2 columns

# Plot all dn and ds for each protein
for (i in 1:length(res)){
  c1 <- res[[i]]
  prot <- names(res)[[i]]
  header <- gsub(".FUBAR.json", "", prot)
  plot(c1$alpha, main = header, col="dodgerblue", pch=19)
  points(c1$beta, main = header, col="orange", pch=19)
}

