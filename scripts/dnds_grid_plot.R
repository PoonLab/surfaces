# Analyse FUBAR results for multiple CDSs in a virus
setwd("/home/laura/Projects/surfaces/scripts_pipeline/second_all_jsons")
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
pdf("/home/laura/Projects/surfaces/evolution_plots/all_dn_ds.pdf",
    width=12, height = 12)
par(mfrow=c(4,4)) # 6 rows, 2 columns
par(mar=c(4,4,2,2)) # 6 rows, 2 columns

# Plot all dn and ds for each protein
for (i in 1:length(res)){
  c1 <- res[[i]]
  prot <- names(res)[[i]]
  header <- gsub(".FUBAR.json", "", prot)
  plot(c1$alpha, main = header, col="dodgerblue", pch=19)
  points(c1$beta, main = header, col="orange", pch=19)
}

dev.off()

###############################################
# Join all dn and ds values
###############################################
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
# breaks <- c(0, 0.1, 0.5, 1, 1.6, 2.5, 5)
breaks <- c(0, 0.1, 0.25, 0.5, 1, 1.6, 2.5, 5, 10, 50)
n <- length(breaks)

pdf("/home/laura/Projects/surfaces/evolution_plots/selected_fingerprints_chik_mumps.pdf", 
    width=12, height = 12)
par(mfrow=c(4,4)) # 6 rows, 2 columns
par(mar=c(4,4,2,2)) # 6 rows, 2 columns
for (i in 1:length(res)){
  name <- strsplit(names(res[i]), "[.]")[[1]][1]
  dnds <- res[[i]]

  x <- findInterval(dnds$alpha, vec=breaks)
  y <- findInterval(dnds$beta, vec=breaks)
  tab <- table(x, y)
  n <- length(breaks)
  
  plot(NA, xlim=c(0,n), ylim=c(0,n), type='n', bty='n',
       xaxt='n', yaxt='n', xlab="Synonymous rates", cex.lab=size.axis,
       ylab="Nonsynonymous rates", 
       # main=names(res)[[i]])
       main=name)
  axis(side=1, at=0:(n-1), labels=c(breaks), cex.axis=size.axis)
  axis(side=2, at=0:(n-1), labels=c(breaks), cex.axis=size.axis, las=1)
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

dev.off()
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
# Functions to annotate dataframe
##################################
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
##################################
# Measure fingerprint similarity
##################################
cosine <- function(x,y) {
  sum(x*y) / sqrt(as.numeric(sum(x*x)) * as.numeric(sum(y*y)))
}

###########################################
# Calculate cosine distance
###########################################
# old breaks
# breaks <- c(0, 0.1, 0.5, 1, 1.6, 2.5, 5)

breaks <- c(0, 0.1, 0.25, 0.5, 1, 1.6, 2.5, 5, 10, 50)

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

#####################################
# Calculate PCA
#####################################
pca <- prcomp(cos.dis)
# plot(pca$x[,1], pca$x[,2])

# Annotate PCA
pca.df <- as.data.frame(pca$x)[, c("PC1", "PC2", "PC3", "PC4")]
pca.df$virus <- unlist(lapply(rownames(pca.df),get_virus))
pca.df$protein <- unlist(lapply(rownames(pca.df),get_protein))

# variation: From Sareh's surfaces/fubar/create_pca.r
pca.var <- pca$sdev^2
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var), 1)


# --------- Annotate PCA with virus family
flaviviridae <- c("zika", "dengue2", "dengue1", "wnv", "YF")
togaviridae <- c("chikv")
bornaviridae <- c("bornediseasevirus")
rhabdoviridae <- c("lyssavirus")
paramyxoviridae <- c("mumps", "measles")


# TO DO: there has to be a better way to do this
pca.df$fam <- ifelse(pca.df$virus %in% flaviviridae, "Flaviviridae",
                     ifelse(pca.df$virus %in% togaviridae, "Togaviridae", 
                            ifelse(pca.df$virus %in% bornaviridae, "Bornaviridae",
                                   ifelse(pca.df$virus %in% rhabdoviridae, "Rhabdoviridae",
                                          ifelse(pca.df$virus %in% paramyxoviridae, "Paramyxoviridae", NA)))))
# --- Annotate surface proteins --- #
# chikv, mumps, lyssavirus, measles, zika, tbe, dengue, yellow fever
surface.proteins <- c("E1_envelope_glycoprotein", "E2_envelope_glycoprotein",
                      "SH_protein", "fusion_protein", "hemagglutinin", "small_hydrophobic_protein", 
                      "glycoprotein", "membrane_glycoprotein_M",
                      "hemagglutinin", "fusion_protein", "hemagglutinin_protein",
                      "envelope_protein_E",
                      "envelope_protein_E",
                      "envelope_protein", "matrix_protein_M", "membrane_glycoprotein_precursor_M"
                      )
pca.df$is.surface <- pca.df$protein %in% surface.proteins

#####################################
# My data Information
#####################################
families <- c("Flaviviridae",
              "Bornaviridae", 
              "Togaviridae",
              "Bornaviridae", 
              "Rhabdoviridae",
              "Paramyxoviridae")

virus <- c("zika","YF", "dengue1", "dengue2","wnv",
           "bornediseasevirus", 
           "chikv",
           "lyssavirus", 
           "measles", "mumps"
           )

#####################################
# Defining colors plots
#####################################
cols <- c("#cf4e4e", "#f48c06", "#f7b801", "#ddb892","#7f5539",
          "#a0b22f",
          "#966b9d", 
          "#ea9ab2",
          "#468faf", "#1e6091"
          # "#adb5bd", "gray", "#539053", "#a8dff7"
          )
names(cols) <- virus

#####################################
# Defining color shapes
#####################################
shapes <- c(21, 22, 23, 24, 25, 3)
names(shapes) <- families

#####################################
# Plot PCA
#####################################

pdf("/home/laura/Projects/surfaces/evolution_plots/pca_disected.pdf", 
    width=8, height = 8)
par(mfrow=c(2,2)) # 6 rows, 2 columns
par(mar=c(3,4,2,2)) # lower, left, upper, right
# Plot PCA colored by virus

plot(pca.df$PC1, pca.df$PC2, las=1,
     cex = 1.4, lwd = 1.5,
     # bg=add.alpha(cols[pca.df$virus], 0.5),
     bg=ifelse(pca.df$is.surface, cols[pca.df$virus], add.alpha(cols[pca.df$virus], 0.2)),
     col = cols[pca.df$virus],
     pch = ifelse(pca.df$is.surface, 22, 21),
     # pch = shapes[pca.df$fam],
     # pch = 19,
     # pch=as.numeric(as.factor(pca.df$fam)),
     xlim = c(-2,4), ylim = c(-4,2),
     xlab = "Principal Component 1", ylab = "Principal Component 2"
     )

legend("bottomleft", legend = names(cols),
       col=cols, pch=19, cex=1.1)

# text(pca.df$PC1[pca.df$is.surface], pca.df$PC2[pca.df$is.surface],
#      labels=pca.df$protein[pca.df$is.surface], pos=4, cex=1)


# --- upper right corner --- #
smaller.cex <- 1.5
smaller.lwd <- 1.8
text.cx <- 0.9

plot(pca.df$PC1, pca.df$PC2, las=1,
     cex = smaller.cex, lwd = smaller.lwd,
     # bg=add.alpha(cols[pca.df$virus], 0.5),
     bg=ifelse(pca.df$is.surface, cols[pca.df$virus], add.alpha(cols[pca.df$virus], 0.2)),
     col = cols[pca.df$virus],
     pch = ifelse(pca.df$is.surface, 22, 21),
     # pch = shapes[pca.df$fam],
     # pch = 19,
     # pch=as.numeric(as.factor(pca.df$fam)),
     xlim = c(1.5,3.5), ylim = c(0,2)
     )

text(pca.df$PC1, pca.df$PC2,
     labels=pca.df$protein, pos=4, cex=text.cx)

# --- lower right corner --- #
plot(pca.df$PC1, pca.df$PC2, las=1,
     cex = smaller.cex, lwd = smaller.lwd,
     # bg=add.alpha(cols[pca.df$virus], 0.5),
     bg=ifelse(pca.df$is.surface, cols[pca.df$virus], add.alpha(cols[pca.df$virus], 0.2)),
     col = cols[pca.df$virus],
     pch = ifelse(pca.df$is.surface, 22, 21),
     # pch = shapes[pca.df$fam],
     # pch = 19,
     # pch=as.numeric(as.factor(pca.df$fam)),
     xlim = c(1.9,2.6), ylim = c(-3.5,-2.8)
     )

text(pca.df$PC1, pca.df$PC2,
     labels=pca.df$protein, pos=4, cex=text.cx)

# --- top left --- #
plot(pca.df$PC1, pca.df$PC2, las=1,
     cex = smaller.cex, lwd = smaller.lwd,
     # bg=add.alpha(cols[pca.df$virus], 0.5),
     bg=ifelse(pca.df$is.surface, cols[pca.df$virus], add.alpha(cols[pca.df$virus], 0.2)),
     col = cols[pca.df$virus],
     pch = ifelse(pca.df$is.surface, 22, 21),
     # pch = shapes[pca.df$fam],
     # pch = 19,
     # pch=as.numeric(as.factor(pca.df$fam)),
     xlim = c(-2,1), ylim = c(-0.5,1)
     )

text(pca.df$PC1, pca.df$PC2,labels=pca.df$protein, pos=4, cex=text.cx)
#      labels=pca.df$protein, pos=4, cex=text.cx)
# 
# text(pca.df$PC1[pca.df$is.surface], pca.df$PC2[pca.df$is.surface],
#      labels=pca.df$protein[pca.df$is.surface], pos=4, cex=text.cx)
# # 
# text(pca.df$PC1, pca.df$PC2,
#      labels=pca.df$protein, pos=4, cex=text.cx)

dev.off()

# Plot variation
barplot(pca.var.per, main="Component contribution", xlab = "PC", las=1,
        ylab="Percent Variation", xlim=c(0,10), col = "dodgerblue4")
legend("topright", legend = unique(pca.df$virus),
       col=as.factor(unique(pca.df$virus)), pch=19, cex=1.1)

#####################################
# UMAP
#####################################
library(umap)

umap.res <- umap(cos.dis)
umap.res <- as.data.frame(umap.res$layout)
colnames(umap.res) <- c("x", "y")
umap.res$virus <- unlist(lapply(rownames(umap.res),get_virus))
umap.res$protein <- unlist(lapply(rownames(umap.res),get_protein))

surface.proteins <- c("E1_envelope_glycoprotein", "E2_envelope_glycoprotein",
                      "SH_protein", "fusion_protein", "hemagglutinin", "small_hydrophobic_protein", 
                      "glycoprotein", "membrane_glycoprotein_M",
                      "hemagglutinin", "fusion_protein",
                      "envelope_protein_E",
                      "envelope_protein_E",
                      "envelope_protein"
)

umap.res$is.surface <- umap.res$protein %in% surface.proteins
pdf("/home/laura/Projects/surfaces/evolution_plots/umap_res.pdf", 
    width=20, height = 10)

par(mfrow=c(1,2)) # 6 rows, 2 columns
# Plot colored by virus
plot(umap.res$x, umap.res$y, col=cols[umap.res$virus], las=1,
     lwd = 1.5, cex =1.2,
     # xlim = c(-4,8),
     xlim=c(-4,2), ylim=c(-2,6),
     # pch = ifelse(umap.res$is.surface, 19, 1)
     pch = shape[]
)

legend("bottomleft", legend = names(cols),
       col=cols, pch=19, cex=0.8)

text(umap.res$x, umap.res$y,
     labels=umap.res$protein, pos=4, cex=0.5)

plot(umap.res$x, umap.res$y, col=as.factor(umap.res$is.surface), 
     # xlim = c(-4,4), ylim = c(-4,4),
     pch=19, las=1 
     # ,main=paste("intervals =", interval)
)
text(umap.res$x[umap.res$is.surface], umap.res$y[umap.res$is.surface],
     labels=umap.res$protein[umap.res$is.surface], pos=4, cex=1)
dev.off()
###################################
# Getting equiprobable breakpoints
###################################
x <- quantile(all.dn, probs = seq(0, 1, length.out = 5))
cut(all.dn, breaks=x)
tiled<-ntile(all.dn, n=5)

###########################################################
# Old attempt: results by changing the number of intervals
##########################################################
# Changing the number of intervals
int <- c(2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 50, 100, 150, 200, 250, 300)
for (j in 1:length(int)) {
  
  # --- Define grid points --- #  
  interval=int[j]
  x <- quantile(all.dn, probs = seq(0, 1, length.out = interval))
  breaks <-x
  
  # --- Do stuff here 
}  

par(mfrow=c(16,3)) #row, column
# par(mfrow=c(2,3)) #row, column
par(mar=c(2, 4, 3, 1))

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


