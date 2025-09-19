require(jsonlite)

breaks <- (50*(0:19)^5)/19^5  # from Murrell et al. 2016


setwd("~/git/surfaces/data/")
#grids <- read_json("fingerprints.json")
mdat <- read.csv("metadata.csv", na.strings="")
mdat$key <- paste(mdat$virus, mdat$protein)

# import alignment stats
astats <- read.csv("align_stats.csv")
# unscramble the rows
idx <- match(mdat$key, paste(astats$virus, astats$protein))
astats <- astats[idx, ]
mdat$ncod <- astats$ncod
mdat$nseq <- astats$nseq
mdat$tree.len <- astats$treelen
mdat$label <- paste(mdat$abbrv, mdat$short)


# generate content for supplementary table
require(knitr)
kable(mdat[order(mdat$family, mdat$abbrv) , 
           c('abbrv', 'protein', 'short', 
             'exposed', #'polymerase', 'protease', 'structural', 
             'ncod', 'nseq', 'tree.len')],
      format='latex', row.names=FALSE, booktabs=TRUE, digits=3, linesep='')

# supplementary figure
pdf("~/papers/surfaces/img/align-stats.pdf", width=9, height=4)
par(mar=c(5,5,1,1), mfrow=c(1,2))
hist(mdat$ncod, main=NA, xlab="Alignment length (codons)")
rug(mdat$ncod)
abline(v=100, lty=2)
hist(mdat$tree.len, main=NA, xlab="Tree length")
rug(mdat$tree.len)
abline(v=c(0.5, 2), lty=2)
dev.off()


# calculate Earth mover's distance
require(transport)
n <- length(breaks)
coords <- expand.grid(1:n, 1:n)

# load `grids` object from get_fingerprints.R
#load("grids.RData")  # L=50, supfig
load("step9_grids.RData")  # L=100

# calculate weighted point pattern for each fingerprint
wpps <- lapply(grids, function(g) {
  wpp(coords, mass=as.numeric(g$grid))
})
virus <- sapply(grids, function(g) g$virus)
protein <- sapply(grids, function(g) gsub(" step[689]", "", g$protein))
names(wpps) <- paste(virus, protein, sep='.')

# these should be the same, but just make sure
idx <- match(paste(virus, protein), mdat$key)
sum(is.na(idx)) == 0  # check
mdatx <- mdat[idx, ]  # expand metadata for replicates

# this takes a minute - about an hour for replicate samples
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
}, mc.cores = 12)


#save(res, file="wasserstein-res-step9.RData")
load("wasserstein-res-step9.RData")  # L=100 codons
#load("wasserstein-res-step6.RData")
#load("wasserstein-res4.RData")  # L=50 codons, sup fig

wmat <- matrix(unlist(res), nrow=n, ncol=n, byrow=T)
# reflect upper triangular portion of matrix
ix <- lower.tri(wmat, diag=FALSE)
wmat[ix] <- t(wmat)[ix]  
rownames(wmat) <- names(wpps)
colnames(wmat) <- names(wpps)  

#write.csv(wmat, file="wdist-sample10.csv", quote=F)
#wmat <- read.csv("wdist-revised.csv", row.names=1)

wdist <- as.dist(wmat)

mds <- cmdscale(wdist, k=2, eig=T)

barplot(head(mds$eig)/sum(mds$eig))
#pal <- ifelse(mdat$protein_classification[mdat$include=="Y"]=="Surface", "red", "cadetblue")
#labels <- gsub("protein", "", names(wpps))
labels <- paste(mdatx$abbrv, mdatx$short)


### STEP 6 ONLY ###
#count.ncod <- read.csv("5_alt/count_ncod.csv")
#idx <- match(names(wpps), count.ncod$key)
#ncod <- count.ncod$ncod[idx]
#pdf("~/papers/surfaces/img/ncodons-raw.pdf", width=5, height=5)
#par(mar=c(0,0,0,0))
#plot(mds$points, cex=sqrt(ncod)/10)
#dev.off()
##text(mds$points[,1], mds$points[,2], labels=names(wpps), cex=0.5)


#pdf("mds.pdf", width=11, height=18)
par(mfrow=c(1,1), mar=c(0,0,0,0))
plot(mds$points[,1:2], type='n')
points(mds$points[,1], mds$points[,2], pch=1, cex=sqrt(mdat$ncod)/10) #, col=pal)
#idx <- order(mdat$ncod)[1:20]
#points(mds[idx, 1], mds[idx,2], cex=2, col='red')
#idx <- which(virus=="IBV")
#text(mds[idx,1], mds[idx,2], labels=labels[idx], cex=0.5) #, labels=mdat$keys[mdat$include=='Y'], cex=0.5) #, col=pal)
text(mds$points[,1], mds$points[,2], labels=labels, cex=0.6,
     col=ifelse(mdatx$exposed, 'red', 'blue')) #, labels=names(wpps), cex=0.5, col=pal)
#dev.off()


# show that replicates are clustered
idx <- as.integer(as.factor(labels))
set.seed(1)
pick <- sample(1:max(idx), 25)
pick.lab <- unique(labels)[pick]

#pdf("~/papers/surfaces/img/samples-cluster.pdf", width=6, height=6)
png("~/papers/surfaces/img/samples-cluster.png", width=6*600, 
    height=6*600, res=600)
par(mar=c(0,0,0,0), mfrow=c(5,5))
for (i in 1:25) {
  plot(mds$points[,1:2], pch=19, col='grey', cex=0.2, main=pick.lab[i], 
       line=-1, xaxt='n', yaxt='n', xlab=NA, ylab=NA, cex.main=1, font.main=1)
  points(mds$points[which(idx==pick[i]), 1], mds$points[which(idx==pick[i]), 2])
}
dev.off()




# calculate centroids for each virus-protein combo
cents <- t(sapply(split(1:nrow(mds$points), mdatx$key), function(i) {
  c(mean(mds$points[i, 1]), mean(mds$points[i, 2]))
}))

idx <- match(row.names(cents), mdat$key)
mdat <- mdat[idx, ]
labels <- paste(mdat$abbrv, mdat$short)


# compare to ncodons-raw.pdf (see above using step 6 data)
pdf("~/papers/surfaces/img/ncodons-adjusted.pdf", width=5, height=5)
par(mar=c(0,0,0,0), mfrow=c(1,1))
plot(cents, cex=sqrt(mdat$ncod)/10)
dev.off()

require(dichromat)

#pdf("~/papers/surfaces/img/mds.pdf", width=10, height=5)
#pdf("~/papers/surfaces/img/mds50.pdf", width=10, height=5)  # sup fig
par(mfrow=c(1,2), mar=c(0,0,0,0))
plot(cents, type='n', bty='n', xaxt='n', yaxt='n')
#text(cents[,1], cents[,2], labels, cex=0.7, col=ifelse(mdat$enveloped, 'red', 'blue'))
text(cents[,1], cents[,2], labels, 
     cex=ifelse(mdat$exposed, 0.7, 0.5), 
     font=ifelse(mdat$exposed, 2, 1),
     col=ifelse(!mdat$exposed, 'cadetblue',
                ifelse(mdat$enveloped, 'firebrick', 'darkorange'))
     )
text(max(cents[,1]), max(cents[,2]), label="Surface-exposed", adj=1, cex=1.2)
abline(v=par('usr')[2])
plot(cents, type='n', yaxt='n', bty='n', xaxt='n', yaxt='n')
text(cents[,1], cents[,2], labels, 
     cex=ifelse(mdat$polymerase, 0.7, 0.5), 
     font=ifelse(mdat$polymerase, 2, 1),
     col=ifelse(!mdat$polymerase, 'cadetblue',
                ifelse(mdat$enveloped, 'firebrick', 'darkorange')), 
)
text(max(cents[,1]), max(cents[,2]), label="Polymerase", adj=1, cex=1.2)
#dev.off()


# calculate silhouette width criterion (SWC) for centroids
#cdmx <- as.matrix(dist(cents))  # centroid distance matrix
cdmx <- matrix(0, nrow=nrow(cents), ncol=nrow(cents))
cdmx[lower.tri(cdmx)] <- dist(cents)
cdmx[upper.tri(cdmx)] <- t(cdmx[lower.tri(cdmx)])

idx <- which(mdat$exposed & mdat$enveloped)
idx2 <- which(mdat$exposed & !mdat$enveloped)
w1 <- cdmx[idx, idx]  # within
swc.a <- apply(w1, 1, sum) / (nrow(w1)-1)
b1 <- cdmx[idx, idx2]  # between
swc.b <- apply (b1, 1, sum) / ncol(b1)

silh <- (swc.b - swc.a) / pmax(swc.a, swc.b)
barplot(silh[order(silh)])
swc <- sum(silh) / length(silh)

w1 <- cdmx[idx2, idx2]  # within
swc.a <- apply(w1, 1, sum) / (nrow(w1)-1)
b1 <- cdmx[idx2, idx]  # between
swc.b <- apply (b1, 1, sum) / ncol(b1)
silh <- (swc.b - swc.a) / pmax(swc.a, swc.b)

############## DEPRECATED ##############

# generate a similar figure that labels enzymes
pdf(file="~/papers/surfaces/img/rdrp.pdf", width=10, height=5)
par(mar=c(0,0,0,0), mfrow=c(1,2))
plot(cents, type='n', yaxt='n')
text(cents[,1], cents[,2], labels, 
     cex=ifelse(mdat$polymerase, 0.7, 0.5), 
     font=ifelse(mdat$polymerase, 2, 1),
     col=ifelse(!mdat$polymerase, 'cadetblue',
                ifelse(mdat$enveloped, 'firebrick', 'darkorange')), 
)
plot(cents, type='n', yaxt='n')
text(cents[,1], cents[,2], labels,
     cex=ifelse(mdat$protease, 0.7, 0.5), 
     font=ifelse(mdat$protease, 2, 1),
     col=ifelse(!mdat$protease, 'cadetblue',
                ifelse(mdat$enveloped, 'firebrick', 'darkorange')), 
)
dev.off()


#points(cents[mdat$enveloped, 1], cents[mdat$enveloped, 2], cex=2)
#idx <- is.element(mdat$virus, c("Foveavirus", "PotatoX", "PotatoY", "Tobacco"))
#text(cents[,1], cents[,2], labels, cex=0.7, col=ifelse(idx, 'red', 'blue'))

plot(cents, cex=sqrt(mdat$ncod)/10)
plot(cents, cex=sqrt(mdat$nseq)/5)
plot(cents, cex=mdat$tree.len)

require(rgl)
mds3 <- cmdscale(wdist, k=3)
cents3 <- t(sapply(split(1:nrow(mds3), mdatx$key), function(i) {
  apply(mds3[i,], 2, mean)
}))

plot3d(cents3, type='n')
text3d(cents3, texts=labels, 
       col=ifelse(!mdat$exposed, 'cadetblue',
                  ifelse(mdat$enveloped, 'firebrick', 'darkorange'))
       )


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
