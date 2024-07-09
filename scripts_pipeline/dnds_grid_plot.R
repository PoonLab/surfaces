# Analyse FUBAR results for multiple CDSs in a virus
setwd("/home/laura/Projects/surfaces/scripts_pipeline/all_jsons")
library(plyr)
require(jsonlite)
library(dplyr)

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

###############################################
# distribution of dN and dS for all proteins
###############################################
par(mfrow=c(5,5)) # 6 rows, 2 columns

# Plot all dn and ds for each protein
for (i in 1:length(res)){
  c1 <- res[[i]]
  prot <- names(res)[[i]]
  header <- gsub(".FUBAR.json", "", prot)
  plot(c1$alpha, main = header, col="dodgerblue", pch=19)
  points(c1$beta, main = header, col="orange", pch=19)
}

par(mfrow=c(1,2)) # 6 rows, 2 columns

# Extract dn ds for all proteins into a single dataframe
all.dn <- unlist(lapply(res, function(d){d$beta}))
all.ds <- unlist(lapply(res, function(d){d$alpha}))
summary(cbind(all.dn, all.ds))
hist(all.dn, breaks=100)
hist(all.ds, breaks = 100)

##################################
# Plot grid for all proteins
##################################
breaks <- c(0, 0.1, 0.5, 1, 1.6, 2.5, 5)
par(mfrow=c(5,5)) # 6 rows, 2 columns
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
###########################################
# cosine distance and PCA for all proteins
###########################################
breaks <- c(0, 0.1, 0.5, 1, 1.6, 2.5, 5)

# --- get cosine distance between all --- #
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
res.names <- unlist(lapply(names(res), function(x) gsub(".FUBAR.json", "", x)))
rownames(cos.dis)<- res.names
colnames(cos.dis)<-res.names

# --- PCA --- #
pca <- prcomp(cos.dis)
plot(pca$x[,1], pca$x[,2])

# Get protein name
get_protein <- function(x){
  splited <- strsplit(x, "_")[[1]]
  protein <- paste(splited[2:length(splited)], collapse = "_")
  return(protein)
}

# Get virus name
get_virus <- function(x){
  splited <- strsplit(x, "_")[[1]]
  virus <- splited[1]
  return(virus)
}

# Annotate PCA
pca.df <- as.data.frame(pca$x)[, c("PC1", "PC2", "PC3", "PC4")]
pca.df$virus <- unlist(lapply(rownames(pca.df),get_virus))
pca.df$protein <- unlist(lapply(rownames(pca.df),get_protein))

# variation: From Sareh's surfaces/fubar/create_pca.r
pca.var <- pca$sdev^2
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var), 1)

# --- Annotate surface proteins --- #
# chikv, mumps, lyssavirus, measles, zika, tbe
surface.proteins <- c("E1_envelope_glycoprotein", "E2_envelope_glycoprotein",
                      "SH_protein", "fusion_protein", "hemagglutinin", "small_hydrophobic_protein", 
                      "glycoprotein",
                      "hemagglutinin", "fusion_protein",
                      "envelope_protein_E",
                      "envelope_protein_E"
                      )

# --- Plot --- #
par(mfrow=c(1,3)) #row, column

# Plot variation
barplot(pca.var.per, main="Component contribution", xlab = "PC", las=1,
        ylab="Percent Variation", xlim=c(0,10), col = "dodgerblue4")

# Plot PCA colored by is.surface
pca.df$is.surface <- pca.df$protein %in% surface.proteins
plot(pca.df$PC1, pca.df$PC2, col=as.factor(pca.df$is.surface), pch=19, las=1)
text(pca.df$PC1[pca.df$is.surface], pca.df$PC2[pca.df$is.surface],
     labels=pca.df$protein[pca.df$is.surface], pos=4, cex=1.1)
legend("topleft", legend = c("not surface", "surface"),
       col=as.factor(unique(pca.df$is.surface)), pch=19)

# Plot PCA colored by virus
plot(pca.df$PC1, pca.df$PC2, col=as.factor(pca.df$virus), pch=19, las=1)
legend("topleft", legend = unique(pca.df$virus),
       col=as.factor(unique(pca.df$virus)), pch=19, cex=1.1)

################################################
# Calculate cosine distance between two clusters
################################################
clust.1 <- res[[1]]
clust.2 <- res[[2]]

fp1 <- get_fingerprint(clust.1, breaks = breaks)
fp2 <- get_fingerprint(clust.2, breaks = breaks)
sim <- cosine(fp1, fp2)
##################################
# Load metadata
##################################
md<-read.csv("feb28measles_protein-clusters-info.csv")
g.n<-table(md$gene.name,md$clusters)
prot.name<-c()
len.clus <-c()
for(i in 1:8){
  m <- which.max(g.n[,i])
  l <- sum(g.n[,6])
  len.clus <- append(len.clus,l)
  prot.name<-append(prot.name,names(m))
}
###################################
# Get tree lengths
###################################
require(ape)
phy.files <- Sys.glob("*.tree")
all.trees <- lapply(phy.files, function(f) read.tree(f))
t.l <- lapply(all.trees, function(t) sum(t$edge.length))
t.t <- lapply(all.trees, function(t) length(t$tip.label))

############################################
# Merge protein information into a dataframe
############################################
prot.d <- data.frame(
              virus = rep("measles", length(phy.files)),
              cluster = c(1:length(phy.files)),
              clus.length = len.clus,
              prot.name = prot.name,
              tree.length = unlist(t.l),
              tree.tips = unlist(t.t),
              is.surface = rep(NA, length(phy.files))
                   ) 

write.csv(prot.d, "measles-surfaces-info.csv", row.names = FALSE)


