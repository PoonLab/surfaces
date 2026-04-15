require(jsonlite)
require(parallel)
require(transport)

# load metadata
setwd("~/git/surfaces/data/")

parse.json <- function(f) {
  # extract information from filename
  tokens <- strsplit(basename(f), "_")[[1]]
  virus <- tokens[1]
  protein <- paste(tokens[2:(length(tokens)-1)], collapse=" ")
  
  # posterior distribution
  json <- read_json(f, simplifyVector = T)
  grid <- matrix(json$grid[,3], nrow=20, byrow=T)
  values <- json$grid[1:20, 2]
  
  list(virus=virus, protein=protein, values=values, grid=grid)
}

# re-ran analysis with full alignments (no normalization of tree lengths)
files <- Sys.glob("iss135/*.json")
grids <- lapply(files, function(f) parse.json(f))

coords <- expand.grid(1:20, 1:20)
wpps <- lapply(grids, function(g) {
  wpp(coords, mass=as.numeric(g$grid))
})

virus <- sapply(grids, function(g) g$virus)
protein <- sapply(grids, function(g) gsub(" step[4689]", "", g$protein))
names(wpps) <- paste(virus, protein, sep='.')
keys <- paste(virus, protein)


# import metadata
mdat <- read.csv("iss135/metadata_iss135.csv", na.strings="")
mdat$key <- paste(mdat$virus, mdat$protein)
idx <- match(keys, mdat$key)
keys[is.na(idx)]
mdat <- mdat[idx, ]

astats <- read.csv("iss135/iss135_align_stats.csv")
astats$key <- paste(astats$virus, astats$protein)
idx <- match(mdat$key, astats$key)

mdat$nseq <- astats$nseq[idx]
mdat$ncod <- astats$ncod[idx]
mdat$treelen <- astats$treelen[idx]

dnds <- read.csv("iss135/dnds.csv")
idx <- match(mdat$key, paste(dnds$virus, dnds$protein)) 
mdat$dnds <- dnds$dnds[idx]
mdat$n.pos <- dnds$n.pos[idx]
mdat$n.neg <- dnds$n.neg[idx]


#### re-do supplementary figure

# original metadata file
mdat0 <- read.csv("metadata.csv", na.strings="")
mdat0$key <- paste(mdat0$virus, mdat0$protein)

# import alignment stats
astats0 <- read.csv("align_stats.csv")
# unscramble the rows
idx <- match(mdat0$key, paste(astats0$virus, astats0$protein))
astats0 <- astats0[idx, ]

# append alignment stats to metadata
mdat0$ncod <- astats0$ncod
mdat0$nseq <- astats0$nseq
mdat0$tree.len <- astats0$treelen
mdat0$label <- paste(mdat0$abbrv, mdat0$short)

corner <- function(px=-0.19, py=1) {
  usr <- par('usr')
  dx <- usr[2] - usr[1]
  dy <- usr[4] - usr[3]
  list(x=px*dx+usr[1], y=py*dy+usr[3])
}



pdf("~/papers/surfaces/img/align-stats.pdf", width=9, height=3)

par(mar=c(5,5,2,1), mfrow=c(1,3))
hist(mdat0$ncod, main=NA, xlab="Alignment length (codons)", border='grey')
text(corner(), label="A", cex=2.5, xpd=NA, font=1)
rug(mdat0$ncod)
abline(v=100, lty=2)

hist(mdat$treelen, main=NA, xlab="Tree length (original)", breaks=20, 
     border='grey')
rug(mdat$treelen)
rect(xl=0.5, xr=2, yb=0, yt=100, col=rgb(1,0,0,0.1), border=NA)
abline(v=c(0.5, 2), col='firebrick')
text(corner(), label="B", cex=2.5, xpd=NA, font=1)

hist(mdat0$tree.len, main=NA, xlab="Tree length (pruned)", border='grey')
rug(mdat0$tree.len)
text(corner(), label="C", cex=2.5, xpd=NA, font=1)
rect(xl=0.5, xr=2, yb=0, yt=100, col=rgb(1,0,0,0.1), border=NA)
abline(v=c(0.5, 2), col='firebrick')

dev.off()



hist(mdat$treelen)
rug(mdat$treelen)

n <- length(wpps)
res <- mclapply(0:(n*n-1), function(k) {
  i <- k %/% n + 1
  j <- k %% n + 1
  if (i < j) {
    wasserstein(wpps[[i]], wpps[[j]], p=2, prob=TRUE)
  } else { 0 }
}, mc.cores = 12)  # this takes a minute

# convert list of distances into matrix
wmat <- matrix(unlist(res), nrow=n, ncol=n, byrow=T)
# reflect upper triangular portion of matrix
ix <- lower.tri(wmat, diag=FALSE)
wmat[ix] <- t(wmat)[ix]  
keys <- gsub("\\.", " ", names(wpps))
rownames(wmat) <- keys
colnames(wmat) <- keys


# do multidimensional scaling
wdist <- as.dist(wmat)
mds <- cmdscale(wdist, k=2, eig=T)

labels <- paste(mdat$abbrv, mdat$short)


pdf("~/papers/surfaces/img/normalization.pdf", width=8, height=4)
par(mar=rep(0.5, 4), mfrow=c(1,2))
o <- order(mdat$ncod, decreasing=T)
plot(mds$points[o,], cex=sqrt(mdat$ncod[o])/10, 
     xaxt='n', yaxt='n', bty='n',
     xlab=NA, ylab=NA, pch=21, col='white', bg='cadetblue',
     main='Alignment length', font.main=1, adj=0, line=-1)

o <- order(mdat$treelen, decreasing=TRUE)
plot(mds$points[o,], cex=0.2+sqrt(mdat$treelen[o])/1.5, xaxt='n', yaxt='n', bty='n',
     xlab=NA, ylab=NA, pch=21, col='white', bg='orange',
     main='Tree length', font.main=1, adj=0, line=-1)
dev.off()


# visualize distribution of tree lengths
plot(mds$points, cex=sqrt(mdat$treelen)/1.5)
text(mds$points, label=labels, cex=0.5)


cor.test(mds$points[,1], mdat$ncod, method='spearman')
cor.test(mds$points[,1], mdat$treelen, method='spearman')
cor.test(mds$points[,2], mdat$treelen, method='spearman')

fit.1 <- glm(mds$points[,1] ~ log(mdat$ncod) + log(mdat$treelen))
z.1 <- residuals(fit.1)
fit.2 <- glm(mds$points[,2] ~ log(mdat$ncod) + log(mdat$treelen))
z.2 <- residuals(fit.2)

# recalculate distances
dmx <- dist(cbind(z.1, z.2))
m2 <- cmdscale(dmx, k=2, eig=F)

plot(m2, cex=sqrt(mdat$treelen), pch=21, bg=rgb(1,0,0,0.3))


par(mar=c(0,0,0,0), mfrow=c(1,1))
plot(m2, col='grey', cex=sqrt(mdat$treelen), xaxt='n', yaxt='n', bty='n')
text(m2, label=labels, cex=0.5, col=ifelse(mdat$exposed, 'red', 'blue'))



pdf("~/papers/surfaces/img/mds-full.pdf", width=10, height=5)
par(mar=c(0,0,0,0), mfrow=c(1,2))
plot(m2, type='n', xaxt='n', yaxt='n', bty='n')
text(m2, label=labels,
     cex=ifelse(mdat$exposed, 0.7, 0.5),
     font=ifelse(mdat$exposed, 2, 1),
     col=ifelse(!mdat$exposed, 'cadetblue', 
                ifelse(mdat$enveloped, 'firebrick', 'darkorange'))
)
text(min(m2[,1]), max(m2[,2]), label="Surface-exposed", adj=0, cex=1.2)
abline(v=par('usr')[2])
plot(m2, type='n', xaxt='n', yaxt='n', bty='n')
text(m2, label=labels, 
     cex=ifelse(mdat$polymerase, 0.7, 0.5), 
     font=ifelse(mdat$polymerase, 2, 1),
     col=ifelse(!mdat$polymerase, 'cadetblue',
                ifelse(mdat$enveloped, 'firebrick', 'darkorange')), 
)
text(min(m2[,1]), max(m2[,2]), label="Polymerase", adj=0, cex=1.2)
dev.off()


# look at association between residualized MDS and gene-wide dN/dS
par(mar=c(5,5,1,1))
plot(m2[,1], mdat$dnds, log='y')
plot(m2[,2], mdat$dnds)
plot(m2[,2], mdat$n.pos/mdat$ncod)

cor.test(m2[,1], mdat$dnds)
cor.test(m2[,2], mdat$dnds)



# display fingerprints at the extremes of either MDS coordinates
left <- c(54, 132, 137) 
right <- c(117, 176, 84)
top <- c(24, 147, 19)
bottom <- c(51, 232, 64)

# change keywords ('top') to make all four components of plot
pdf("~/papers/surfaces/img/grids-top.pdf", width=9, height=3)
par(mfrow=c(1,3), mar=c(1,1,1,1))
#pal <- colorRampPalette(c('white', 'midnightblue'))(18)
pal <- colorRampPalette(c('white', 'red4'))(18)
pal <- c('white', pal[3:18])
for (i in top) {
  image(grids[[i]]$grid, col=pal, xaxt='n', yaxt='n')
  abline(a=0, b=1)
  abline(v=0.5, lty=2)
  abline(h=0.5, lty=2)
  text(x=0.03, y=0.95, adj=0, label=paste(mdat$abbrv[i], mdat$short[i]), cex=2)
}
dev.off()




idx <- which(!mdat$enveloped)
plot(m2, type='n')
text(m2[idx,], label=labels[idx], cex=0.6, col=ifelse(mdat$exposed[idx], 'red', 'blue'))


### SAVE OUR WORK
#save(m2, labels, mdat, dmx, mds, wmat, file="iss135_final3.RData")
load("iss135_final3.RData")



d.ncod <- as.matrix(dist(log(mdat$ncod)))
x <- d.ncod[upper.tri(d.ncod)]
y <- wmat[upper.tri(wmat)]

require(MASS)

png("~/papers/surfaces/img/logcodon.png", width=5*300, height=5*300, res=300)
par(mar=c(5,5,1,1))
plot(x, y, bty='n', xlab="log(codons)", ylab="Wasserstein distance",
     cex=0.5, col='black', lwd=2.5, xaxt='n')
axis(side=1, )
points(x, y, col='cadetblue3', pch=19, cex=0.5)
contour(kde2d(x, y), add=TRUE)
dev.off()

cor.test(x, y)

# what if we don't log-transform the alignment lengths?
d.ncod <- as.matrix(dist(mdat$ncod))
x <- d.ncod[upper.tri(d.ncod)]
cor.test(x, y)

#write.csv(as.matrix(dmx), file="iss135-dmx2.csv", quote=F)

d.tlen <- as.matrix(dist(log(mdat$treelen)))
x <- d.tlen[upper.tri(d.tlen)]
plot(x, y, pch=19, cex=0.3)



#######################################

# show the corrective effects of residualization or downsampling



load("L100.RData")  # overwrites `wmat` and `mdatx`
wdist <- as.dist(wmat)
mds <- cmdscale(wdist, k=2, eig=T)

# calculate centroids for each virus-protein combo
cents <- t(sapply(
  split(1:nrow(mds$points), factor(mdatx$key, levels=unique(mdatx$key))), 
  function(i) {
    c(mean(mds$points[i, 1]), mean(mds$points[i, 2]))
  }
))

# reconstitute one-per-row metadata
mdat2 <- mdatx[!duplicated(mdatx), ]


# m2 should be loaded from above
pdf("~/papers/surfaces/img/corrected-mds.pdf", width=7, height=7)
par(mfrow=c(2,2), mar=c(1,2,2,1))

idx <- order(mdat$ncod, decreasing=TRUE)
plot(m2[idx,], cex=sqrt(mdat$ncod[idx])/10, 
     main="Alignment length", adj=0, font.main=1, 
     bty='n', xaxt='n', yaxt='n', xlab=NA, ylab=NA,
     pch=21, col='white', bg='cadetblue')

title(ylab="Residualized", line=1, cex.lab=1.2)

idx <- order(mdat$treelen, decreasing=TRUE)
plot(m2[idx,], cex=0.2+sqrt(mdat$treelen[idx])/1.5,
     main="Tree length", adj=0, font.main=1,
     pch=21, col='white', bg='orange',
     bty='n', xaxt='n', yaxt='n', xlab=NA, ylab=NA)

par(mar=c(2,2,1,1))

idx <- order(mdat2$ncod, decreasing=TRUE)
plot(cents[idx,], cex=sqrt(mdat2$ncod[idx])/10, 
     bty='n', xaxt='n', yaxt='n', xlab=NA, ylab=NA, 
     pch=21, col='white', bg='cadetblue')

abline(h=par('usr')[4]*1.1, xpd=NA, col='grey')
title(ylab="Downsampled", line=1, cex.lab=1.2)

idx <- order(mdat2$tree.len, decreasing=TRUE)
plot(cents, cex=0.2+sqrt(mdat2$tree.len)/1.5, 
     pch=21, col='white', bg='orange',
     bty='n', xaxt='n', yaxt='n', xlab=NA, ylab=NA)

dev.off()


#######################

# Compare distance matrices

# average distances among replicates
load("L50.RData")  # provides `wmat` and `mdatx`

idxs <- sapply(split(1:nrow(mdatx), mdatx$key), function(x) {
  if (length(x)==1) {
    rep(x, 10)
  } else { x }
})

# reconstitute distance matrices
res <- lapply(1:10, function(z) {
  idx <- idxs[z, ]
  wmat[idx, idx]
})

wmat.avg <- Reduce("+", res) / length(res)
write.csv(as.matrix(wmat.avg), file="iss135-avgmat-L50.csv", quote=F)


###############################

m1 <- as.matrix(read.csv("iss135-dmx2.csv", row.names=1))
m2 <- as.matrix(read.csv("iss135-avgmat-L100.csv", row.names=1))

common <- intersect(colnames(m1), colnames(m2))
idx <- match(common, colnames(m1))
m1c <- m1[idx, idx]
idx <- match(common, colnames(m2))
m2c <- m2[idx, idx]

require(ade4)
res <- mantel.test(m1c, m2c, nperm=1e4, graph=TRUE)
res <- mantel.randtest(as.dist(m1c), as.dist(m2c), nrepet=99999)


x <- m1c[upper.tri(m1c)]
y <- m2c[upper.tri(m2c)]
plot(x, y, pch=19, cex=0.5, col=rgb(0,0,0,0.2))
cor.test(x, y)


#############################

load("iss135_final3.RData")

require(vegan)
adonis2(dmx ~ exposed, data=mdat)
adonis2(dmx ~ exposed * enveloped, data=mdat, permutations=99999, by='terms')



pdf("~/papers/surfaces/img/mds-full.pdf", width=8, height=4)

par(mfrow=c(1,2), mar=c(0,0,0,0))
plot(m2, type='n', xaxt='n', yaxt='n', bty='n', xlab=NA, ylab=NA,
     main=' Surface-exposed', font.main=1, adj=0, line=-1)
abline(v=par('usr')[2])
points(m2, cex=ifelse(mdat$exposed, 0, 0.3), pch=19, col='grey')
text(m2, label=labels,
     cex=ifelse(mdat$exposed, 0.6, 0.01),
     font=ifelse(mdat$enveloped, 3, 2),
     col=ifelse(!mdat$enveloped, 'cadetblue4', 'darkorange2')
)

text(x=-1.3, y=-1.05, label="Non-enveloped", cex=1,
     col='cadetblue4', font=2)
text(x=-0.1, y=-1.05, label="Enveloped", cex=1, col='darkorange2', font=3)

plot(m2, type='n', xaxt='n', yaxt='n', bty='n', xlab=NA, ylab=NA,
     main=' Non-exposed', font.main=1, adj=0, line=-1)
points(m2, cex=ifelse(mdat$exposed, 0.3, 0), pch=19, col='grey')
text(m2, label=labels,
     cex=ifelse(!mdat$exposed, 0.6, 0.01),
     font=ifelse(mdat$enveloped, 3, 2),
     col=ifelse(!mdat$enveloped, 'cadetblue4', 'darkorange2')
)
dev.off()


fam <- c("Flaviviridae", "Orthomyxoviridae", "Paramyxoviridae", 
         "Picornaviridae", "Retroviridae", "Togaviridae")
pal <- gg.rainbow(n=6, l=60)

pdf("~/papers/surfaces/img/by-family.pdf", width=6, height=4)
par(mfrow=c(2,3), mar=c(0,0,0,0))
for (i in 1:6) {
  f <- fam[i]
  plot(m2, type='n', xaxt='n', yaxt='n', bty='n', xlab=NA, ylab=NA)
  points(m2, pch=ifelse(mdat$exposed, 19, 21), 
         cex=ifelse(mdat$family==f, 1.5, 0.3),
         col=ifelse(mdat$family==f, pal[i], 'grey'))
  title(main=paste(" ", f), font.main=1, adj=0, line=-1)
}
dev.off()


adonis2(dmx~family + exposed + enveloped, data=mdat, by='terms', permutations=1e5)

idx <- which(mdat$family=='Retroviridae')  #'Picornaviridae')
temp <- as.dist(as.matrix(dmx)[idx, idx])
adonis2(temp ~ exposed, data=mdat[idx,], permutations=1e5, by='terms')

restrict <- how(blocks=mdat$family, nperm=1e4-1)
adonis2(dmx ~ exposed*family, data=mdat, by='terms', permutations=restrict)



