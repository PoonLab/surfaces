# how to partition dN/dS space
breaks = c(0, 0.1, 0.25, 0.5, 1, 1.6, 2.5, 5, 10, 50)

fingerprint <- function(dnds) {
  x <- findInterval(dnds$alpha, vec=breaks)
  y <- findInterval(dnds$beta, vec=breaks)
  xx <- factor(x, levels=1:9)
  yy <- factor(y, levels=1:9)
  tab <- table(xx, yy)
  tab/sum(tab)  # normalize
}

# from Hugo
cosine <- function(x, y) {
  1 - sum(x * y) / sqrt(sum(x * x) * sum(y * y))  
}

setwd("~/git/surfaces")
fubar <- read.csv("data/collated-fubar.csv")
fubar <- fubar[fubar$step==5, ]

by.gene <- split(fubar, f=list(fubar$protein, fubar$virus), drop=TRUE)
prints <- lapply(by.gene, fingerprint)

# cosine distances
n <- length(prints)
csd <- matrix(0, nrow=n, ncol=n)
for (i in 1:(n-1)) {
  for (j in (i+1):n) {
    csd[i,j] <- csd[j,i] <- cosine(prints[[i]], prints[[j]])
  }
}


# visualize
mds <- cmdscale(csd)
par(mfrow=c(1,1), mar=c(0,0,0,0))
plot(mds, type='n')
text(mds[,1], mds[,2], label=names(by.gene), cex=0.4)

hc <- hclust(as.dist(csd))
plot(hc)


# repeat with Hugo's file
setwd("~/git/surfaces")
info <- read.csv("data/surfaces_info.csv")
info <- info[grepl("step5", info$filename), ]
by.gene <- split(info, f=list(info$filename), drop=TRUE)
prints <- lapply(by.gene, fingerprint)
n <- length(prints)

# cosine distance matrix
csd <- matrix(0, nrow=n, ncol=n)
for (i in 1:(n-1)) {
  for (j in (i+1):n) {
    csd[i,j] <- csd[j,i] <- cosine(prints[[i]], prints[[j]])
  }
}
mds <- cmdscale(csd)
par(mfrow=c(1,1), mar=c(0,0,0,0))
plot(mds, type='n')
text(mds[,1], mds[,2], label=names(by.gene), cex=0.4)
is.surface <- sapply(split(info$protein_classification, info$filename), 
                     function(x) x[1] == "Surface")
points(mds[,1], mds[,2], pch=19, col=ifelse(is.surface, 'blue', 'red'))


# annotate by tree length

# (venv) art@aziraphale:~/git/surfaces/data/5_alt$ for f in *.fasta; 
# do fasttree -nt -quote $f > ${f%.fasta}.nwk; done
require(ape)
df <- info[info$pos==1, ]
df$nwkfile <- paste(df$virus_name, '_', df$gene, '_step5.nwk', sep='')
setwd("~/git/surfaces/data/5_alt")
df$tree.len <- sapply(df$nwkfile, function(f) {
  phy <- read.tree(f)
  sum(phy$edge.length)
})
# had to manually rename a lot of files.....

z <- log(df$tree.len)
z <- (z-min(z)) / (max(z) - min(z))

par(mar=c(0,0,0,0))
plot(mds, type='n')
points(mds[,1], mds[,2], pch=19, col=rgb(z, 0.5, 1-z), cex=z+0.2)


par(mar=c(5,5,1,1))
boxplot(split(df$tree.len, df$virus_name), las=2, cex.axis=0.8,
        ylab='Tree length')


#rownames(csd) <- paste(df$virus_name, df$gene)
rownames(csd) <- names(by.gene)
hc <- hclust(as.dist(csd), method='ward.D2')

require(dendextend)
dend <- as.dendrogram(hc)
dend <- set(dend, "labels_cex", 0.4)

png("~/Desktop/dendogram.png", width=5*300, height=15*300, res=300)
par(mar=c(0,0,0,5))
plot(dend, horiz=T)
dev.off()



# let's try earth mover's distance instead, reusing code from tragula
require(transport)

# let bins be defined by the midpoint between log(breaks)
mid.pts <- log10(diff(breaks)/2 + breaks[1:(length(breaks)-1)])
coords <- expand.grid(mid.pts, mid.pts)

# calculate weighted point pattern for each fingerprint
wpps <- lapply(prints, function(fp) {
  wpp(coords, mass=as.numeric(fp))
})

# this takes a minute
n <- length(wpps)
res <- lapply(0:(n*n-1), function(k) {
  i <- k %/% n + 1
  j <- k %% n + 1
  if (i < j) {
    wasserstein(wpps[[i]], wpps[[j]], p=2, prob=TRUE)
  } else {
    0
  }
})

wdist <- matrix(unlist(res), nrow=n, ncol=n, byrow=T)
# reflect upper triangular portion of matrix
ix <- lower.tri(wdist, diag=FALSE)
wdist[ix] <- t(wdist)[ix]  
rownames(wdist) <- names(wpps)
colnames(wdist) <- names(wpps)  
wdist <- as.dist(wdist)

par(mar=c(0, 0, 0, 0))
mds <- cmdscale(wdist, k=2)
labels <- gsub("_step5\\.fubar", "", names(wpps))

plot(mds[,1:2], type='n')
#text(mds[,1], mds[,2], labels=labels, cex=0.5)
text(mds[,1], mds[,2], cex=0.5,
     labels=sapply(labels, function(x) strsplit(x, "_")[[1]][1]))

# can re-use vector from cosine analysis
points(mds[,1], mds[,2], pch=19, 
       cex=ifelse(is.surface, 1, 0.5),
       col=ifelse(is.surface, 'dodgerblue', 'firebrick'))


pdf("hclust.pdf", width=15, height=5)
hc <- hclust(wdist)
par(mar=rep(0, 4))
plot(hc, cex=0.5, main=NA)
dev.off()


# I'm concerned about using one set of breaks for both alpha and beta
par(mar=c(5,5,1,1))
hist(log10(info$alpha), freq=F, col=rgb(0,0,1,0.5), main=NA, breaks=100, border=NA)
hist(log10(info$beta), col=rgb(1,0,0,0.5), add=T, freq=F, breaks=100, border=NA)

# try different sets
nb <- 7
a.breaks <- 10^quantile(log10(info$alpha), prob=seq(0, 1, length.out=nb))
b.breaks <- 10^quantile(log10(info$beta), prob=seq(0, 1, length.out=nb))
abline(v=log10(a.breaks), col='blue', lwd=0.5)

# manual tuning
hist(log10(info$alpha), col=rgb(0,0,1,0.5), main=NA, breaks=100, border=NA)
a.breaks <- 10^c(-3, -0.2, 0, 0.2, 0.4, 0.5, 0.7, 0.85, 2)
abline(v=log10(a.breaks), col='blue')

hist(log10(info$beta), col=rgb(1,0,0,0.5), main=NA, breaks=100, border=NA)
b.breaks <- 10^c(-5, -1.9, -1.5, -1.2, -0.9, -0.5, -0.2, 0.2, 2)
abline(v=log10(b.breaks), col='red')

fprint <- function(dnds) {
  x <- findInterval(dnds$alpha, vec=a.breaks)
  y <- findInterval(dnds$beta, vec=b.breaks)
  xx <- factor(x, levels=1:(length(a.breaks)-1))
  yy <- factor(y, levels=1:(length(b.breaks)-1))
  tab <- table(xx, yy)
  tab/sum(tab)  # normalize
}

prints <- lapply(by.gene, fprint)

a.mids <- log10(diff(a.breaks)/2 + a.breaks[1:(length(a.breaks)-1)])
b.mids <- log10(diff(b.breaks)/2 + b.breaks[1:(length(b.breaks)-1)])
coords <- expand.grid(a.mids, b.mids)

wpps <- lapply(prints, function(fp) {
  wpp(coords, mass=as.numeric(fp))
})
n <- length(wpps)
res <- mclapply(0:(n*n-1), function(k) {
  i <- k %/% n + 1
  j <- k %% n + 1
  if (i < j) {
    wasserstein(wpps[[i]], wpps[[j]], p=2, prob=TRUE)
  } else {
    0
  }
}, mc.cores=4)
wdist <- matrix(unlist(res), nrow=n, ncol=n, byrow=T)
ix <- lower.tri(wdist, diag=FALSE)
wdist[ix] <- t(wdist)[ix]  
rownames(wdist) <- names(wpps)
colnames(wdist) <- names(wpps)  
wdist <- as.dist(wdist)

mds <- cmdscale(wdist, k=2)
labels <- gsub("_step5\\.fubar", "", names(wpps))
par(mar=c(0, 0, 0, 0))
plot(mds[,1:2], type='n')
points(mds[,1], mds[,2], pch=19, 
       cex=ifelse(is.surface, 1, 0.5),
       col=ifelse(is.surface, 'dodgerblue', 'firebrick'))
text(mds[,1], mds[,2], labels=labels, cex=0.5)


pdf("hclust2.pdf", width=15, height=5)
hc <- hclust(wdist)
par(mar=rep(0, 4))
plot(hc, labels=labels, cex=0.5, main=NA)
dev.off()

