require(jsonlite)

breaks <- (50*(0:19)^5)/19^5  # from Murrell et al. 2016


setwd("~/git/surfaces/data/")
#grids <- read_json("fingerprints.json")
mdat <- read.csv("metadata.csv", na.strings="")

# calculate Earth mover's distance
require(transport)
n <- length(breaks)
coords <- expand.grid(1:n, 1:n)

# calculate weighted point pattern for each fingerprint
wpps <- lapply(grids, function(g) {
  wpp(coords, mass=as.numeric(g$grid))
})
virus <- sapply(grids, function(g) g$virus)
protein <- sapply(grids, function(g) g$protein)
names(wpps) <- paste(virus, protein, sep='.')

# these should be the same, but just make sure
idx <- match(paste(virus, protein), paste(mdat$virus, mdat$protein))
mdat <- mdat[idx, ]

# this takes a minute
require(parallel)
n <- length(wpps)
res <- mclapply(0:(n*n-1), function(k) {
  i <- k %/% n + 1
  j <- k %% n + 1
  if (i < j) {
    wasserstein(wpps[[i]], wpps[[j]], p=2, prob=TRUE)
  } else {
    0
  }
}, mc.cores = 10)

wmat <- matrix(unlist(res), nrow=n, ncol=n, byrow=T)
# reflect upper triangular portion of matrix
ix <- lower.tri(wmat, diag=FALSE)
wmat[ix] <- t(wmat)[ix]  
rownames(wmat) <- names(wpps)
colnames(wmat) <- names(wpps)  
#write.csv(wmat, file="wdist-revised.csv", quote=F)
write.csv(wmat, file="wdist-clean.csv", quote=F)

#wmat <- read.csv("wdist-revised.csv", row.names=1)

wdist <- as.dist(wmat)

mds <- cmdscale(wdist, k=2)
#pal <- ifelse(mdat$protein_classification[mdat$include=="Y"]=="Surface", "red", "cadetblue")
#labels <- gsub("protein", "", names(wpps))
labels <- paste(mdat$abbrv, mdat$short)

pdf("mds.pdf", width=11, height=18)
par(mfrow=c(1,1), mar=c(0,0,0,0))
plot(mds[,1:2], type='n')
#points(mds[,1], mds[,2], pch=19) #, col=pal)
#idx <- which(virus=="IBV")
#text(mds[idx,1], mds[idx,2], labels=labels[idx], cex=0.5) #, labels=mdat$keys[mdat$include=='Y'], cex=0.5) #, col=pal)
text(mds[,1], mds[,2], labels=labels, cex=0.5,
     col=ifelse(mdat$exposed, 'red', 'blue')) #, labels=names(wpps), cex=0.5, col=pal)
dev.off()


require(rgl)
mds3 <- cmdscale(wdist, k=3)
plot3d(mds3)
text3d(mds3, texts=labels, col=ifelse(mdat$exposed, 'red', 'blue'))


require(uwot)
u <- uwot::umap(wdist, n_components=2)

pdf("umap.pdf", width=8, height=8)
par(mar=c(0,0,0,0))
plot(u, type='n')
#labels <- gsub("nonstructural|protein", "", row.names(u))
text(u[,1], u[,2], labels=labels, cex=0.7,
     col=ifelse(mdat$exposed, 'red', 'blue'))#, col=pal)
dev.off()

hc <- hclust(wdist)
pdf("hclust.pdf", width=15, height=5)
par(mar=c(0,0,0,0))
plot(hc, cex=0.5, main=NA, labels=names(wpp), col=pal)
dev.off()


# generate a graph
par(mar=c(5,5,1,1))
hist(wmat[upper.tri(wmat)], breaks=100, xlim=c(0, 1))
rug(wmat[upper.tri(wmat)])
table(wmat[upper.tri(wmat)] < 1.0)

#cutoff <- max(apply(wmat, 1, function(x) min(x[x>0])))
cutoff <- 1.
n <- nrow(wmat)
nm <- gsub("protein|non.*structural", "", names(wpps)) #names(wmat))
k <- 3
adj <- matrix(0, nrow=n, ncol=n, dimnames=list(nm, nm))
for (i in 1:n) {
  row <- as.numeric(wmat[i,])
  ranks <- order(row)
  nn <- ranks[-(ranks==i)][1:k]
  nn <- nn[row[nn] <= cutoff]
  adj[i, nn] <- 1
}

require(igraph)
g <- graph_from_adjacency_matrix(adj, mode="max")
pdf("graph.pdf", width=8, height=8)
plot(g, #layout=layout_with_gem, 
     vertex.label.cex=0.8, vertex.size=2, 
     vertex.color=pal, vertex.border=NA)
dev.off()

write.csv(adj, file="adjmat.csv", quote=FALSE)
