## Requires wmat and mdatx from earthmover.R
setwd("~/git/surfaces/data")
load("L100.RData")  # provides `wmat` and `mdatx`
#load("L50.RData")

# make sure labels match up
keys <- gsub("\\.", " ", rownames(wmat))
stopifnot(all(keys==mdatx$keys))

# do multidimensional scaling
wdist <- as.dist(wmat)
mds <- cmdscale(wdist, k=2, eig=T)
all(rownames(mds$points) == mdatx$key) 
labels <- paste(mdatx$abbrv, mdatx$short)

###### USED FOR STEP 6 DATA ONLY ######
## MDS is skewed by alignment length (number of codons)
# count.ncod <- read.csv("5_alt/count_ncod.csv")
# idx <- match(names(wpps), count.ncod$key)
# ncod <- count.ncod$ncod[idx]
# pdf("~/papers/surfaces/img/ncodons-raw.pdf", width=5, height=5)
# par(mar=c(0,0,0,0))
# plot(mds$points, cex=sqrt(ncod)/10)
# dev.off()
#### can be replaced with:
#par(mfrow=c(1,1), mar=c(0,0,0,0))
#plot(mds$points, cex=sqrt(mdatx$ncod)/10)
#######################################


# show that replicates are clustered
idx <- as.integer(as.factor(labels))
set.seed(1)
pick <- sample(1:max(idx), 25)
pick.lab <- unique(labels)[pick]

##pdf("~/papers/surfaces/img/samples-cluster.pdf", width=6, height=6)
#png("~/papers/surfaces/img/samples-cluster.png", width=6*600, 
#    height=6*600, res=600)
par(mar=c(0,0,0,0), mfrow=c(5,5))
for (i in 1:25) {
  plot(mds$points[,1:2], pch=19, col='grey', cex=0.2, xaxt='n', yaxt='n', 
       xlab=NA, ylab=NA)
  title(main=pick.lab[i], line=-1, font.main=1, cex.main=1)
  points(mds$points[which(idx==pick[i]), 1], mds$points[which(idx==pick[i]), 2])
}
#dev.off()


# Look at just one replicate (this only works for L=100)
idx <- seq(4, nrow(wmat), 10)
wm2 <- wmat[idx, idx]
mds2 <- cmdscale(as.dist(wm2))
mdat <- mdatx[idx,]
par(mar=c(0,0,0,0))
plot(mds2, type='n')
text(mds2, label=mdat$label, cex=ifelse(mdat$exposed, 0.7, 0.5), 
     font=ifelse(mdat$exposed, 2, 1),
     col=ifelse(!mdat$exposed, 'cadetblue', 
                ifelse(mdat$enveloped, 'firebrick', 'darkorange'))
)


# calculate centroids for each virus-protein combo
cents <- t(sapply(
  split(1:nrow(mds$points), factor(mdatx$key, levels=unique(mdatx$key))), 
  function(i) {
    c(mean(mds$points[i, 1]), mean(mds$points[i, 2]))
    }
  ))

# reconstitute one-per-row metadata
mdat <- mdatx[!duplicated(mdatx), ]
# make sure the metadata still matches up with MDS data
stopifnot(all(rownames(cents) == mdat$key))
# which(rownames(cents) != mdat$key)
labels <- paste(mdat$abbrv, mdat$short)

# compare to ncodons-raw.pdf (see above using step 6 data)
#pdf("~/papers/surfaces/img/ncodons-adjusted.pdf", width=5, height=5)
par(mar=c(0,0,0,0), mfrow=c(1,1))
plot(cents, cex=sqrt(mdat$ncod)/10)
#dev.off()

#require(dichromat)

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


## Look at flaviviruses
# flavi <- unique(mdat$abbrv[mdat$family=='Flaviviridae'])
# pdf("~/papers/surfaces/img/flavi.pdf", width=5, height=7.5)
# par(mfrow=c(3,2), mar=c(0,0,0,0))
# for (ab in flavi) {
#   idx <- which(mdat$abbrv==ab)
#   plot(cents, type='n', xaxt='n', yaxt='n')
#   text(cents[idx, 1], cents[idx,2], labels=labels[idx], cex=1,
#        col=ifelse(grepl("^C", mdat$short[idx]), "red", "grey50"),
#        font=ifelse(grepl("^C", mdat$short[idx]), 2, 1))  
# }
# dev.off()


####################

require(kernlab)

labels <- ifelse(mdat$enveloped[mdat$exposed], 1, -1)

# use the first replicate for L=100 samples (step 9)
idx <- seq(4, nrow(wmat), 10)
wm1 <- wmat[idx, idx]  # distances
gamma <- 1 # 1
xm1 <- as.matrix(gamma * exp(-wm1[mdat$exposed, mdat$exposed]^2))
#max.val <- max(xm1)
#xm1 <- (max.val-xm1)/(2*max.val)  # convert to similarity on (0,1)

# k-fold cross-validation
pos <- which(labels>0)
pos <- sample(pos, length(pos), replace=F)  # permutation
pos <- split(pos, ceiling(seq_along(pos)/10))
neg <- which(labels<0)
neg <- sample(neg, length(neg), replace=F)
neg <- split(neg, ceiling(seq_along(neg)/7))
for (k in 1:3) {
  test <- c(neg[[k]], pos[[k]])
  train <- seq_along(labels)[c(-neg[[k]], -pos[[k]])]
  fit <- ksvm(x=as.kernelMatrix(xm1[train, train]), 
              y=labels[train], kernel='matrix', type='C-svc')  
  pred <- predict(fit, newdata=as.kernelMatrix(xm1[test, train]))
  print(table(labels[test], pred))
}
# this doesn't work.......


# leave-one-out
for (j in 1:10) {
  idx <- seq(j, nrow(wmat), 10)
  wm1 <- wmat[idx, idx]  # distances
  xm1 <- as.matrix(exp(-wm1[mdat$exposed, mdat$exposed]^2))  # RBF kernel
  
  res <- matrix(0, nrow=2, ncol=2)
  for (i in 1:length(labels)) {
    fit <- ksvm(x=as.kernelMatrix(xm1[-i, -i]),
                y=labels[-i], kernel='matrix', type='C-svc')
    pred <- predict(fit, newdata=as.kernelMatrix(xm1[i, -i, drop=FALSE]))
    j <- max(labels[i]+1, 1)
    k <- max(pred+1, 1)
    res[j,k] <- res[j,k] + 1
  }
  print(res)  
}



##############################################################

# Try using SVM on centroids of MDS
mds <- cmdscale(wdist, k=10, eig=T)
cents <- t(sapply(split(1:nrow(mds$points), mdatx$key), function(i) {
  apply(mds$points[i, ], 2, mean)
  #c(mean(mds$points[i, 1]), mean(mds$points[i, 2]))
}))

temp <- cents[mdat$exposed, ]
n <- nrow(temp)
labels <- ifelse(mdat$enveloped[mdat$exposed], 1, -1)

# generate balanced k-fold test/train subsets
#train <- sample(1:n, round(n*0.67))
k <- 2
pos <- sample(which(labels>0))
p.grps <- split(pos, rep(1:k, length.out=length(pos)))
neg <- sample(which(labels<0))
n.grps <- split(neg, rep(1:k, length.out=length(neg)))

for (i in 1:k) {
  train <- c(p.grps[[i]], n.grps[[i]])
  fit <- ksvm(x=temp[train, ], y=labels[train], kernel='rbfdot', type="C-svc")
  pred <- predict(fit, temp[-train, ])
  tab <- table(labels[-train], pred)
  print(tab)
}

# try leave-one-out validation instead
cnfsn <- matrix(0, nrow=2, ncol=2)
for (i in 1:n) {
  fit <- ksvm(x=temp[-i, ], y=labels[-i], kernel='rbfdot', type="C-svc")
  pred <- predict(fit, temp[i, , drop=FALSE])
  ci <- (labels[i]+1)/2 + 1
  cj <- (pred+1)/2 + 1
  cnfsn[ci, cj] <- cnfsn[ci, cj] + 1
}

2*cnfsn[2,2] / (2*cnfsn[2,2] + cnfsn[1,2] + cnfsn[2,1])  # F1 score

# I'm not sure how valid it is to use SVM on the results of MDS....


############################################################

# try hierachical clustering directly on Wasserstein distances
hc <- hclust(wdist, method='ward.D2')
plot(hc, labels=FALSE)
x <- cutree(hc, k=2)  # two clusters
table(x, mdatx$exposed, mdatx$enveloped)

x <- cutree(hc, k=3)  # three clusters
table(x, mdatx$exposed, mdatx$enveloped)
