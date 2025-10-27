############################# SETTINGS #############################

## uncomment one of the following:
ncodons <- 50  # set for step 8 data (L=50 codons)
#ncodons <- 100  # set for step 9 data (L=100 codons)
# ncodons <- Inf  # set for step 6 data (no down-sampling)

outfile <- "L50.RData"

####################################################################

# load metadata
setwd("~/git/surfaces/data/")
mdat <- read.csv("metadata.csv", na.strings="")
mdat$key <- paste(mdat$virus, mdat$protein)

# import alignment stats
astats <- read.csv("align_stats.csv")
# unscramble the rows
idx <- match(mdat$key, paste(astats$virus, astats$protein))
astats <- astats[idx, ]

# append alignment stats to metadata
mdat$ncod <- astats$ncod
mdat$nseq <- astats$nseq
mdat$tree.len <- astats$treelen
mdat$label <- paste(mdat$abbrv, mdat$short)


## generate content for supplementary table
#require(knitr)
# kable(mdat[order(mdat$family, mdat$abbrv) , 
#            c('abbrv', 'protein', 'short', 
#              'exposed', #'polymerase', 'protease', 'structural', 
#              'ncod', 'nseq', 'tree.len')],
#       format='latex', row.names=FALSE, booktabs=TRUE, digits=3, linesep='')


## supplementary figure
## pdf("~/papers/surfaces/img/align-stats.pdf", width=9, height=4)
# par(mar=c(5,5,1,1), mfrow=c(1,2))
# hist(mdat$ncod, main=NA, xlab="Alignment length (codons)")
# rug(mdat$ncod)
# abline(v=100, lty=2)
# hist(mdat$tree.len, main=NA, xlab="Tree length")
# rug(mdat$tree.len)
# abline(v=c(0.5, 2), lty=2)
## dev.off()


# load `grids` object from get_fingerprints.R
if (ncodons == 50) {
  load("grids.RData")  # L=50, supfig - not some of these are not replicated 10x
} else if (ncodons==100) {
  load("step9_grids.RData")  # L=100  
} else {
  load("step6_grids.RData")  # original flavour
}


# calculate weighted point pattern for each fingerprint
require(transport)
#breaks <- (50*(0:19)^5)/19^5  # from Murrell et al. 2016
coords <- expand.grid(1:20, 1:20)
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


## Uncomment to calculate Earth movers distance de novo
# require(parallel)
# res <- mclapply(0:(n*n-1), function(k) {
#   i <- k %/% n + 1
#   j <- k %% n + 1
#   if (i < j) {
#     wasserstein(wpps[[i]], wpps[[j]], p=2, prob=TRUE)
#   } else {
#     0
#   }  # this takes a minute / an hour for replicate samples
# }, mc.cores = 12)


## Used to store results
# save(res, file="wasserstein-res-step9.RData")

## Load previous calculations (`res` matrix)
if (ncodons==100) {
  load("wasserstein-res-step9.RData")  # L=100 codons  
} else if (ncodons==50) {
  load("wasserstein-res4.RData")  # L=50 codons, sup fig  
} else {
  load("wasserstein-res-step6.RData")  # without down-sampling codons 
}


# convert list of distances into matrix
n <- length(wpps)
stopifnot(length(res) == n*n)
wmat <- matrix(unlist(res), nrow=n, ncol=n, byrow=T)
# reflect upper triangular portion of matrix
ix <- lower.tri(wmat, diag=FALSE)
wmat[ix] <- t(wmat)[ix]  
keys <- gsub("\\.", " ", names(wpps))
rownames(wmat) <- keys
colnames(wmat) <- keys

stopifnot(nrow(wmat) == nrow(mdatx))
stopifnot(all(keys == mdatx$key))

# `mdatx` and `wmat` will now have the same rows
# proceed to mds.R or svm.R
save(wmat, mdatx, file=outfile)
