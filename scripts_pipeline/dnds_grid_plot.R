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
# Calculate cosine distance
##################################
library(magrittr)
library(Hmisc)

# Function cut2: https://rdrr.io/cran/Hmisc/src/R/cut2.s
cuts.alpha <- cut2(clust.2$alpha, g = 6, onlycuts = T) # determine cuts based on df1
cuts.beta <-  cut2(clust.2$beta, g = 6, onlycuts = T) # determine cuts based on df1

cut2(df1$x, cuts = cuts) %>% table
cut2(df2$x, cuts = cuts) %>% table*2 # multiplied by two for better comparison

cut()

breaks <- c(0, 0.1, 0.5, 1, 1.6, 2.5, 5)
c.a <- cuts.alpha
c.b <- cuts.beta

sel.matrix <- function(x, ca, cb){
  ds <- findInterval(x$alpha, vec=c.a) # cuts based on alpha dist
  dn <- findInterval(x$beta, vec=c.b) # cuts based on beta dist
  # Fingerprint: frequency of dn and ds in a grid
  tab <- table(ds, dn)
  return(as.matrix(tab))
}

clust.1 <- res[[1]]
sel.1 <- sel.matrix(clust.1, cuts.alpha, cuts.beta)
clust.2 <- res[[2]]
sel.2 <- sel.matrix(clust.2, cuts.alpha, cuts.beta)


cosineSimilarity <- function(x,y) {
  sim <- x %*% t(y)/(sqrt(rowSums(x^2) %*% t(rowSums(y^2))))
  return(sim)
}

cos.sim <- cosineSimilarity(sel.1, sel.2)
cos.dis <- 1 - cos.sim
pca <- prcomp(cos.dis)

plot(pca$x[,1], pca$x[,2])
