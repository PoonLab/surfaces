# Analyse FUBAR results for multiple CDSs in a virus
setwd("/home/laura/Projects/surfaces_data")
library(plyr)
require(jsonlite)

######################################################
# Function for transparency on plots
######################################################
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

##################################
# Import JSON files 
##################################

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

##################################
# Load metadata
##################################
md<-read.csv("feb28measles_protein-clusters-info.csv")
g.n<-table(md$gene.name,md$clusters)
prot.name<-c()
for(i in 1:8){
  m <- which.max(g.n[,i])
  prot.name<-append(prot.name,names(m))
}

##################################
# Plot grid for all proteins
##################################

breaks <- c(0, 0.1, 0.5, 1, 1.6, 2.5, 5)
par(mfrow=c(2,4)) # 6 rows, 2 columns
for (i in 1:length(res)){
  # name <- strsplit(names(res[i]), "[.]")[[1]][1]
  dnds <- res[[i]]

  x <- findInterval(dnds$alpha, vec=breaks)
  y <- findInterval(dnds$beta, vec=breaks)
  tab <- table(x, y)
  
  par(mar=c(5,5,1,1))
  plot(NA, xlim=c(0,7), ylim=c(0,7), type='n', bty='n',
       xaxt='n', yaxt='n', xlab="Synonymous rates",
       ylab="Nonsynonymous rates", main=prot.name[i])
  axis(side=1, at=0:7, labels=c(breaks, 10), cex.axis=0.8)
  axis(side=2, at=0:7, labels=c(breaks, 10), cex.axis=0.8, las=1)
  for (i in 1:dim(tab)[1]) {
    for (j in 1:dim(tab)[2]) {
      z <- (1 - tab[i,j]/max(tab))^3
      if (z<1) {
        z <- min(0.95, z)
      }
      rect(xl=i, xr=i+1, yb=j-1, yt=j, 
           col=rgb(z,z,z), border=NA)
    }
  }
  abline(a=0, b=1, lty=2)
}

##################################
# Create fingerprint
##################################

get_fingerprint <- function(x, breaks){
  n <- length(breaks)
  ds <- factor(findInterval(x$alpha, vec=breaks), levels=1:n)
  dn <- factor(findInterval(x$beta, vec=breaks), levels=1:n)
  fp <- table(dn, ds) # Fingerprint
  return(fp)
}

##################################
# Measure fingerprint similarity
##################################

cosine <- function(x,y) {
  sum(x*y) / sqrt(as.numeric(sum(x*x)) * as.numeric(sum(y*y)))
}

################################################
# Calculate cosine distance between two clusters
################################################

breaks <- c(0, 0.1, 0.5, 1, 1.6, 2.5, 5)

clust.1 <- res[[1]]
clust.2 <- res[[2]]

fp1 <- get_fingerprint(clust.1, breaks = breaks)
fp2 <- get_fingerprint(clust.2, breaks = breaks)
sim <- cosine(fp1, fp2)

###################################
# cosine distance for all proteins
###################################
prot.sim <- matrix(NA, length(res), length(res))

for (i in 1:length(res)){
  c1 <- res[[i]]
  fp1 <- get_fingerprint(c1, breaks = breaks)
  
  for (j in 1:length(res)) {
    c2 <- res[[j]]
    fp2 <- get_fingerprint(c2, breaks = breaks)
    sim <- cosine(fp1, fp2)
    prot.sim[i,j] <- sim
  }
}

cos.dis <- 1 - prot.sim
pca <- prcomp(cos.dis)
plot(pca$x[,1], pca$x[,2])

###################################
# Get tree lengths
###################################
require(ape)
phy.files <- Sys.glob("*.tree")
all.trees <- lapply(phy.files, function(f) read.tree(f))
t.l <- lapply(all.trees, function(t) sum(t$edge.length))

############################################
# Merge protein information into a dataframe
############################################
prot.d <- data.frame(
              virus = rep("measles", length(phy.files)),
              cluster = c(1:length(phy.files)),
              prot.name = prot.name,
              tree.length = unlist(t.l),
              is.surface = rep(NA, length(phy.files))
                   ) 

write.csv(prot.d, "measles-surfaces-info.csv", row.names = FALSE)

