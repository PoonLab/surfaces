# This script reproduces the evolutionary fingerprinting analysis in 
# Munoz-Baena, Castelan Sanchez et al.

require(jsonlite)
require(parallel)
require(transport)
require(vegan)
require(ggfree)
require(MASS)

dest.dir <- "~/Downloads/"  # change for your filesystem
out.dir <- "~/Desktop/"     # where to write image files
n.cores <- parallel::detectCores() - 1
use.logoffset <- TRUE

# CUSTOM FUNCTIONS

# locate upper left corner of plot region
corner <- function(px=-0.19, py=1) {
  usr <- par('usr')
  dx <- usr[2] - usr[1]
  dy <- usr[4] - usr[3]
  list(x=px*dx+usr[1], y=py*dy+usr[3])
}

# Parse JSON FUBAR output files
parse.json <- function(f) {
  tokens <- strsplit(basename(f), "_")[[1]]
  virus <- tokens[1]
  protein <- paste(tokens[2:(length(tokens)-1)], collapse=" ")
  json <- read_json(f, simplifyVector = TRUE)
  grid <- matrix(json$grid[, 3], nrow=20, byrow=TRUE)
  values <- json$grid[1:20, 2]
  list(virus=virus, protein=protein, values=values, grid=grid)
}


################################
##         LOAD DATA          ##
################################

# Obtain FUBAR JSON files from Zenodo (>300MB)
download.file("https://zenodo.org/records/20753078/files/step4_json.tar.gz?download=1",
  destfile=file.path(dest.dir, "step4_json.tar.gz"), method='wget', quiet=T)
untar(file.path(dest.dir, "step4_json.tar.gz"))
files <- Sys.glob("*_step4.fubar.json")
grids <- lapply(files, parse.json)


# Construct weighted point patterns for fingerprints
if (use.logoffset) {
  log.offset <- log(grids[[1]]$values + 0.05)
  coords <- expand.grid(log.offset, log.offset)[c(2,1)]  
} else {
  # NOTE: some plots generated below assume log-offset coordinates
  coords <- expand.grid(1:20, 1:20)[c(2,1)]  
}

wpps <- lapply(grids, function(g) {
  wpp(coords, mass=as.numeric(g$grid))
})

# Extract virus and protein from filenames
virus <- sapply(grids, function(g) g$virus)
protein <- sapply(grids, function(g) gsub(" step[4689]", "", g$protein))
names(wpps) <- paste(virus, protein, sep='.')
keys <- paste(virus, protein)


# match metadata to JSON data
mdat <- read.csv("https://raw.githubusercontent.com/PoonLab/surfaces/refs/heads/main/data/metadata.csv", 
                 na.strings="")
mdat$key <- paste(mdat$virus, mdat$protein)
idx <- match(keys, mdat$key)
mdat <- mdat[idx, ]

# load alignment statistics
astats <- read.csv("https://raw.githubusercontent.com/PoonLab/surfaces/refs/heads/main/data/align_stats.csv")

# append alignment statistics to metadata
idx <- match(mdat$key, paste(astats$virus, astats$protein))
mdat$nseq <- astats$nseq[idx]
mdat$ncod <- astats$ncod[idx]
mdat$treelen <- astats$treelen[idx]

# append gene-wide selection results
dnds <- read.csv("https://raw.githubusercontent.com/PoonLab/surfaces/refs/heads/main/data/dnds.csv")
idx <- match(mdat$key, paste(dnds$virus, dnds$protein)) 
mdat$dnds <- dnds$dnds[idx]
mdat$n.pos <- dnds$n.pos[idx]
mdat$n.neg <- dnds$n.neg[idx]



###################################################
##  Figure 1 was drawn manually                  ##
##  Use scripts/plot-dnds.R to generate Figure 2 ##
###################################################


#######################################
#  Generate components for Figure 3  ##
#######################################

panels <- list(left = c(54, 132, 137), 
               right = c(117, 176, 84),
               top = c(24, 147, 19),
               bottom = c(51, 232, 64))
for (pname in names(panels)) {
  if (pname %in% c("top", "right")) {
    pal <- colorRampPalette(c('white', 'midnightblue'))(18)
  } else {
    pal <- colorRampPalette(c('white', 'red4'))(18)
  }
  pal <- c('white', pal[3:18])
  
  pdf(file.path(out.dir, paste0("grids-", pname, ".pdf")), width=9, height=3)
  par(mfrow=c(1, 3), mar=c(1, 1, 1, 1))
  for (i in panels[[pname]]) {
    image(grids[[i]]$grid, col=pal, xaxt='n', yaxt='n')
    abline(a=0, b=1)
    abline(v=0.5, lty=2)
    abline(h=0.5, lty=2)
    text(x=0.03, y=0.95, adj=0, label=paste(mdat$abbrv[i], mdat$short[i]), cex=2)
  }
  dev.off()
}



#############################################
##     Calculate Wasserstein distances     ##
#############################################

n <- length(wpps)
res <- mclapply(0:(n*n-1), function(k) {
  i <- k %/% n + 1
  j <- k %% n + 1
  if (i < j) {
    wasserstein(wpps[[i]], wpps[[j]], p=2, prob=TRUE)
  } else { 0 }
}, mc.cores = n.cores)  # this takes a minute on an AMD Ryzen 9 7950X

# convert list of distances into matrix
wmat <- matrix(unlist(res), nrow=n, ncol=n, byrow=T)
ix <- lower.tri(wmat, diag=FALSE)  # reflect upper portion
wmat[ix] <- t(wmat)[ix]  
keys <- gsub("\\.", " ", names(wpps))
rownames(wmat) <- keys
colnames(wmat) <- keys

# do multidimensional scaling
wdist <- as.dist(wmat)
mds <- cmdscale(wdist, k=2, eig=T)


##################################################
# Fingerprints are confounded by alignment characteristics

# correlation between Wasserstein distances and log-transformed alignment lengths
d.ncod <- as.matrix(dist(log(mdat$ncod)))
x <- d.ncod[upper.tri(d.ncod)]
y <- wmat[upper.tri(wmat)]
cor.test(x, y)

# Supplementary Figure S3
png(file.path(out.dir, "logcodon.png"), width=5*300, height=5*300, res=300)
par(mar=c(5,5,1,1))
plot(x, y, bty='n', xlab="Abs. diff. log(codons)", ylab="Wasserstein distance",
     cex=0.5, col='black', lwd=2.5, xaxt='n')
axis(side=1, )
points(x, y, col='cadetblue3', pch=19, cex=0.5)
contour(kde2d(x, y), add=TRUE)
dev.off()


# Supplementary Figure S5
pdf(file.path(out.dir, "normalization.pdf"), width=8, height=4)
par(mar=rep(0.5, 4), mfrow=c(1,2))
o <- order(mdat$ncod, decreasing=T)
plot(mds$points[o,], cex=sqrt(mdat$ncod[o])/10, 
     xaxt='n', yaxt='n', bty='n', xpd=NA,
     xlab=NA, ylab=NA, pch=21, col='white', bg='cadetblue')
title(main='Alignment length', font.main=1, adj=0, line=-1)
points(x=c(-3, -2.5, -1.9, -1.2), y=rep(-2, 4), 
       cex=sqrt(c(50, 100, 500, 1000))/10, pch=21, bg='white', col='black')
text(x=c(-3, -2.5, -1.9, -1.2), y=c(-1.9, -1.88, -1.85, -1.81), 
     labels=c(50, 100, 500, 1000), cex=0.6)
o <- order(mdat$treelen, decreasing=TRUE)
plot(mds$points[o,], cex=0.2+sqrt(mdat$treelen[o])/1.5, xaxt='n', yaxt='n', bty='n',
     xlab=NA, ylab=NA, pch=21, col='white', bg='orange', xpd=NA)
title(main='Tree length', font.main=1, adj=0, line=-1)
points(x=c(-3, -2.6, -2.2, -1.7, -1.1, -0.3), y=rep(-2, 6), xpd=NA,
       cex=0.2+sqrt(c(1, 2, 5, 10, 20, 50))/1.5, pch=21, bg='white', col='black')
text(x=c(-3, -2.6, -2.2, -1.7, -1.1, -0.3), y=c(-1.9, -1.88, -1.86, -2, -2, -2), 
     labels=c(1, 2, 5, 10, 20, 50), cex=0.6)
dev.off()

# technical confounders are correlated with MDS coordinates
cor.test(mds$points[,1], mdat$ncod, method='spearman')
cor.test(mds$points[,2], mdat$ncod, method='spearman')
cor.test(mds$points[,1], mdat$treelen, method='spearman')
cor.test(mds$points[,2], mdat$treelen, method='spearman')

# residualization
fit.1 <- glm(mds$points[,1] ~ log(mdat$ncod) + log(mdat$treelen))
z.1 <- residuals(fit.1)
fit.2 <- glm(mds$points[,2] ~ log(mdat$ncod) + log(mdat$treelen))
z.2 <- residuals(fit.2)

# recalculate distances
dmx <- dist(cbind(z.1, z.2))
m2 <- cmdscale(dmx, k=2, eig=F)


##########################
#     MAKE FIGURE 4      #
##########################

labels <- paste(mdat$abbrv, mdat$short)
pdf(file.path(out.dir, "Figure4.pdf"), width=7, height=3.5)
par(mfrow=c(1,2), mar=c(0,0,0,0))
plot(m2, type='n', xaxt='n', yaxt='n', bty='n', xlab=NA, ylab=NA,
     main=' Non-enveloped', font.main=1, adj=0, line=-1)
abline(v=par('usr')[2])
points(m2, cex=ifelse(mdat$enveloped, 0.5, 0), pch=19, col='grey')
text(m2, label=labels, cex=ifelse(!mdat$enveloped, 0.6, 0.01),
     font=ifelse(!mdat$exposed, 3, 2),
     col=ifelse(mdat$exposed, 'darkorange3', 'orange'))
text(x=-1.5, y=-1.05, label="Exposed", cex=0.8, col='darkorange3', font=2)
text(x=-0.4, y=-1.05, label="Non-exposed", cex=0.8, col='orange', font=3)
rect(xl=-1.95, xr=0.25, yb=-1.115, yt=-0.975, lwd=0.5)
plot(m2, type='n', xaxt='n', yaxt='n', bty='n', xlab=NA, ylab=NA,
     main=' Enveloped', font.main=1, adj=0, line=-1)
points(m2, cex=ifelse(mdat$enveloped, 0, 0.5), pch=19, col='grey')
idx <- order(mdat$exposed)
text(m2[idx, ], label=labels[idx], cex=ifelse(mdat$enveloped[idx], 0.6, 0.01),
     font=ifelse(!mdat$exposed[idx], 3, 2),
     col=ifelse(mdat$exposed[idx], 'steelblue4', 'slategray3'))
text(x=-1.5, y=-1.05, label="Exposed", cex=0.8, col='steelblue4', font=2)
text(x=-0.4, y=-1.05, label="Non-exposed", cex=0.8, col='slategray3', font=3)
rect(xl=-1.95, xr=0.25, yb=-1.115, yt=-0.975, lwd=0.5)
dev.off()


# PERMANOVA analysis directly on Wasserstein distance matrix
adonis2(wdist ~ log(ncod) + log(treelen) + exposed * enveloped, 
        data = mdat, permutations = 99999, by = "terms")
adonis2(wdist ~ log(ncod) + log(treelen) + polymerase * enveloped, 
        data = mdat, permutations = 99999, by = "terms")
adonis2(wdist ~ log(ncod) + log(treelen) + family + exposed * enveloped, 
        data = mdat, permutations = 99999, by = "terms")

restrict <- how(blocks=mdat$family, nperm=1e4-1)
adonis2(wdist ~ log(ncod) + log(treelen) + exposed*family, data=mdat, by='terms', 
        permutations=restrict)


#############################################
##  Family-wise PERMANOVA tests (Table 2)  ##
#############################################

# one versus rest
n.fam <- length(unique(mdat$family))
one.vs <- data.frame(family=character(n.fam), R2=numeric(n.fam), P=numeric(n.fam))
for (i in 1:n.fam) {
  fam <- unique(mdat$family)[i]
  pmn <- adonis2(wmat ~ log(ncod) + log(treelen) + I(family==fam), data=mdat, 
                 by='terms', permutations=9999) 
  one.vs[i, ] <- list(family=fam, R2=pmn$R2[3], P=pmn$`Pr(>F)`[3])
}
one.vs$p.adj <- p.adjust(one.vs$P, method="BH")
one.vs <- one.vs[order(one.vs$family), ]
one.vs$R2 <- round(100*one.vs$R2, 2)

# within (exposed)
within <- data.frame(family=character(n.fam), R2=numeric(n.fam), P=numeric(n.fam))
for (i in 1:n.fam) {
  fam <- unique(mdat$family)[i]
  idx <- which(mdat$family==fam)
  wd <- as.dist(as.matrix(wdist)[idx, idx])
  md <- mdat[idx, ]
  pmn <- adonis2(wd ~ log(ncod) + log(treelen) + exposed, data=md, 
                 permutations=9999, by='terms')
  within[i, ] <- list(family=fam, R2=pmn$R2[3], P=pmn$`Pr(>F)`[3])
}
within$p.adj <- p.adjust(within$P, method="BH")
within <- within[order(within$family), ]
within$R2 <- round(100*within$R2, 2)



# Supplementary Figure S10
fam <- c("Flaviviridae", "Orthomyxoviridae", "Paramyxoviridae", 
         "Picornaviridae", "Retroviridae", "Togaviridae")
pal <- gg.rainbow(n=6, l=60)
pdf(file.path(out.dir, "by-family.pdf"), width=6, height=4)
par(mfrow=c(2,3), mar=c(0,0,0,0))
for (i in 1:6) {
  f <- fam[i]
  plot(m2, type='n', xaxt='n', yaxt='n', bty='n', xlab=NA, ylab=NA)
  points(m2[mdat$family!=f, ], col='grey', pch=19, cex=0.3)
  
  temp <- m2[mdat$family==f, ]
  is.exposed <- mdat$exposed[mdat$family==f]
  points(temp[is.exposed, ], pch=21, bg=pal[i], col='white', cex=2)
  points(temp[!is.exposed, ], pch=21, col=pal[i], cex=1.5, lwd=1.5)
  title(main=paste(" ", f), font.main=1, adj=0, line=-1)
}
dev.off()



# load down-sampled data
download.file("https://zenodo.org/records/20801064/files/L100_logoffset.RData?download=1",
              destfile="L100_logoffset.RData")
load("L100_logoffset.RData")  # CAUTION: THIS OVERWRITES `wmat`!
wdist.1 <- as.dist(wmat)
mds.1 <- cmdscale(wdist.1, k=2, eig=T)
cents <- t(sapply(
  split(1:nrow(mds.1$points), factor(mdatx$key, levels=unique(mdatx$key))), 
  function(i) {  # calculate centroids for each virus-protein combo
    c(mean(mds.1$points[i, 1]), mean(mds.1$points[i, 2]))
  }
))
mdat2 <- mdatx[!duplicated(mdatx), ]  # reconstitute one-per-row metadata
astat0 <- read.csv("https://raw.githubusercontent.com/PoonLab/surfaces/refs/heads/main/data/step7_align_stats.csv")
idx <- match(mdat2$key, paste(astat0$virus, astat0$protein))
mdat2$ncod <- astat0$ncod[idx]
mdat2$tree.len <- astat0$treelen[idx]
wmat <- as.matrix(wdist)  # restore previous `wmat`


## Supplementary Figure S7
pdf(file.path(out.dir, "corrected-mds.pdf"), width=7, height=7)
par(mfrow=c(2,2), mar=c(1,2,2,1))
idx <- order(mdat$ncod, decreasing=TRUE)
plot(m2[idx,], cex=sqrt(mdat$ncod[idx])/10, 
     main="Alignment length", adj=0, font.main=1, 
     bty='n', xaxt='n', yaxt='n', xlab=NA, ylab=NA,
     pch=21, col='white', bg='cadetblue')
title(ylab="Residualized", line=1, cex.lab=1.2)
points(x=c(-0.97, -0.85, -0.7, -0.5), y=rep(0.45, 4), xpd=NA, 
       cex=sqrt(c(50, 100, 500, 1000))/10, pch=21, bg='white', col='black')
text(x=c(-0.97, -0.85, -0.7, -0.5), y=c(.425, .42, .405, .395), 
     labels=c(50, 100, 500, 1000), cex=0.6, xpd=NA)
idx <- order(mdat$treelen, decreasing=TRUE)
plot(m2[idx,], cex=0.2+sqrt(mdat$treelen[idx])/1.5,
     main="Tree length", adj=0, font.main=1,
     pch=21, col='white', bg='orange',
     bty='n', xaxt='n', yaxt='n', xlab=NA, ylab=NA)
points(x=c(-0.93, -0.83, -0.72, -0.61, -0.48, -0.3), y=rep(0.45, 6), xpd=NA,
       cex=0.2+sqrt(c(1, 2, 5, 10, 20, 50))/1.5, pch=21, bg='white', col='black')
text(x=c(-0.93, -0.83, -0.72, -0.61, -0.48, -0.3), 
     y=c(0.415, 0.41, 0.405, 0.45, 0.45, 0.45), 
     labels=c(1, 2, 5, 10, 20, 50), cex=0.6, xpd=NA)
par(mar=c(2,2,1,1))
idx <- order(mdat2$ncod, decreasing=TRUE)
plot(cents[idx,], cex=sqrt(mdat2$ncod[idx])/10, 
     bty='n', xaxt='n', yaxt='n', xlab=NA, ylab=NA, 
     pch=21, col='white', bg='cadetblue')
points(x=c(0.25, 0.32, 0.4, 0.5), y=rep(0.3, 4), xpd=NA, 
       cex=sqrt(c(50, 100, 500, 1000))/10, pch=21, bg='white', col='black')
text(x=c(0.25, 0.32, 0.4, 0.5), y=c(0.282, 0.278, 0.27, 0.26), 
     labels=c(50, 100, 500, 1000), cex=0.6, xpd=NA)
abline(h=par('usr')[4]*1.1, xpd=NA, col='grey')
title(ylab="Downsampled", line=1, cex.lab=1.2)
idx <- order(mdat2$tree.len, decreasing=TRUE)
plot(cents, cex=0.2+sqrt(mdat2$tree.len)/1.5, 
     pch=21, col='white', bg='orange',
     bty='n', xaxt='n', yaxt='n', xlab=NA, ylab=NA)
dev.off()



# association with mode of transmission
mdat$trans.mode <- ifelse(
  mdat$vector, 'vector', ifelse(
    mdat$chronic, 'STBBI', ifelse(
      mdat$fecaloral, 'fecal-oral', ifelse(
        mdat$contact, 'contact', 'respiratory')
    )))
mdat$trans.mode[mdat$abbrv=='BDV'] <- NA
mdat$respiratory <- (!is.na(mdat$trans.mode) & mdat$trans.mode=='respiratory')

# fingerprints are significantly separated by mode of transmission
adonis2(wmat ~ log(ncod) + log(treelen) + trans.mode, data=mdat,
        by='terms', na.action=na.omit, permutations=99999)
adonis2(wmat ~ log(ncod) + log(treelen) + trans.mode*exposed, data=mdat,
        by='terms', na.action=na.omit, permutations=99999)


######################################################
##      Supplementary Table S2
load("L100_logoffset.RData")  # CAUTION: THIS OVERWRITES `wmat`!
idxs <- sapply(split(1:nrow(mdatx), mdatx$key), function(x) {
  if (length(x)==1) { rep(x, 10) } else { x }
})
# reconstitute distance matrices
res <- lapply(1:10, function(z) {
  idx <- idxs[z, ]
  wmat[idx, idx]
})
wmat.100 <- Reduce("+", res) / length(res)  # averaged (centroids)
mdat.100 <- mdatx[idxs[1,], ]
wmat <- as.matrix(wdist)

adonis2(wmat.100 ~ exposed*enveloped, data=mdat.100, permutations=9999, by='terms')
adonis2(wmat.100 ~ family, data=mdat.100, permutations=9999, by='terms')


mds.100 <- cmdscale(wmat.100)
idx <- order(mdat.100$exposed)

##  Supplementary Figure S8
pdf(file.path(out.dir, "mds100.pdf"), width=7, height=3.5)
par(mfrow=c(1,2), mar=c(0,0,0,0))
plot(mds.100, type='n', xaxt='n', yaxt='n', bty='n', xlab=NA, ylab=NA)
abline(v=par('usr')[2])
title(main=' Non-enveloped', font.main=1, adj=0, line=-1)
title(main="100 codons ", font.main=1, adj=1, line=-1)
points(mds.100, cex=ifelse(mdat.100$enveloped, 0.3, 0), pch=19, col='grey')
text(mds.100[idx,], label=mdat.100$label[idx],
     cex=ifelse(!mdat.100$enveloped[idx], 0.6, 0.01),
     font=ifelse(!mdat.100$exposed[idx], 3, 2),
     col=ifelse(mdat.100$exposed[idx], 'darkorange3', 'orange'))
rect(xl=0.29, xr=1.05, yb=-0.64, yt=-0.48, xpd=NA, col='white', lwd=0.5)
text(x=1, y=-0.53, label="Exposed", cex=0.9, col='darkorange3', font=2, 
     xpd=NA, adj=1)
text(x=1, y=-0.6, label="Non-exposed", cex=0.9, col='orange', font=3,
     xpd=NA, adj=1)
plot(mds.100, type='n', xaxt='n', yaxt='n', bty='n', xlab=NA, ylab=NA)
title(main='  Enveloped', font.main=1, adj=0, line=-1)
title(main="100 codons ", font.main=1, adj=1, line=-1)
points(mds.100, cex=ifelse(mdat.100$enveloped, 0, 0.3), pch=19, col='grey')
text(mds.100[idx,], label=mdat.100$label[idx],
     cex=ifelse(mdat.100$enveloped[idx], 0.6, 0.01),
     font=ifelse(!mdat.100$exposed[idx], 3, 2),
     col=ifelse(mdat.100$exposed[idx], 'steelblue4', 'slategray3'))
rect(xl=0.29, xr=1.05, yb=-0.64, yt=-0.48, xpd=NA, col='white', lwd=0.5)
text(x=1, y=-0.53, label="Exposed", cex=0.9, col='steelblue4', font=2, 
     xpd=NA, adj=1)
text(x=1, y=-0.6, label="Non-exposed", cex=0.9, col='slategray3', font=3,
     xpd=NA, adj=1)
dev.off()

###################################################

download.file("https://zenodo.org/records/20801064/files/L50_logoffset.RData?download=1",
              destfile="L50_logoffset.RData")
load("L50_logoffset.RData")
idxs <- sapply(split(1:nrow(mdatx), mdatx$key), function(x) {
  if (length(x)==1) { rep(x, 10) } else { x }
})
res <- lapply(1:10, function(z) {idx <- idxs[z, ]; wmat[idx, idx] })
wmat.50 <- Reduce("+", res) / length(res) 
mdat.50 <- mdatx[idxs[1,], ]
wmat <- as.matrix(wdist)  # restore original matrix

adonis2(wmat.50 ~ exposed*enveloped, data=mdat.50, permutations=9999, by='terms')
adonis2(wmat.50 ~ family, data=mdat.50, permutations=9999, by='terms')


## Supplementary Figure S9
mds.50 <- cmdscale(wmat.50)
idx <- order(mdat.50$exposed)

pdf(file.path(out.dir, "mds50.pdf"), width=7, height=3.5)
par(mfrow=c(1,2), mar=c(0,0,0,0))
plot(mds.50, type='n', xaxt='n', yaxt='n', bty='n', xlab=NA, ylab=NA)
title(main=' Non-enveloped', font.main=1, adj=0, line=-1)
title(main="50 codons ", font.main=1, adj=1, line=-1)
points(mds.50, cex=ifelse(mdat.50$enveloped, 0.3, 0), pch=19, col='grey')
text(mds.50[idx, ], label=mdat.50$label[idx],
     cex=ifelse(!mdat.50$enveloped[idx], 0.6, 0.01),
     font=ifelse(!mdat.50$exposed[idx], 3, 2),
     col=ifelse(mdat.50$exposed[idx], 'darkorange3', 'orange'))
abline(v=par('usr')[2], lwd=0.5)
rect(xl=-1.03, xr=-0.51, yb=-0.42, yt=-0.33, xpd=NA, col='white', lwd=0.5)
text(x=-1, y=-0.36, label="Exposed", cex=0.8, col='darkorange3', 
     font=2, xpd=NA, adj=0)
text(x=-1, y=-0.4, label="Non-exposed", cex=0.8, col='orange', 
     font=3, xpd=NA, adj=0)
plot(mds.50, type='n', xaxt='n', yaxt='n', bty='n', xlab=NA, ylab=NA)
title(main='  Enveloped', font.main=1, adj=0, line=-1)
title(main="50 codons ", font.main=1, adj=1, line=-1)
points(mds.50, cex=ifelse(mdat.50$enveloped, 0, 0.3), pch=19, col='grey')
text(mds.50[idx, ], label=mdat.50$label[idx],
     cex=ifelse(mdat.50$enveloped[idx], 0.6, 0.01),
     font=ifelse(!mdat.50$exposed[idx], 3, 2),
     col=ifelse(mdat.50$exposed[idx], 'steelblue4', 'slategray3'))
rect(xl=0.27, xr=0.79, yb=-0.42, yt=-0.33, xpd=NA, col='white', lwd=0.5)
text(x=0.75, y=-0.36, label="Exposed", cex=0.8, col='steelblue4', 
     font=2, xpd=NA, adj=1)
text(x=0.75, y=-0.4, label="Non-exposed", cex=0.8, col='slategray3', 
     font=3, xpd=NA, adj=1)
dev.off()




######################################################
##  Supplementary Table S3

# 1. No plant proteins exposed = all FALSE
# 2. coat proteins are exposed = T, F, T, F, F, T, F, T
# 3. VSRs are exposed = F, T, F, T, T, F, T, T
# 4. both coat and VSR exposed = all TRUE
# 5. no plant viruses = all NA

plants <- mdat
opts <- list(
  "No plant proteins exposed" = rep(F, 8),
  "Coat proteins are exposed" = c(T, F, T, F, F, T, F, T),
  "VSRs are exposed" = c(F, T, F, T, T, F, T, T),
  "Both coat and VSR exposed" = rep(T, 8),
  "No plant viruses" = rep(NA, 8))
for (i in 1:length(opts)) {
  print(names(opts)[i])
  opt <- opts[[i]]
  plants$exposed[plants$abbrv=="PVX" & plants$short=="CP"] <- opt[1]
  plants$exposed[plants$abbrv=="PVX" & plants$short=="TGB1"] <- opt[2]  # P25
  plants$exposed[plants$abbrv=="PVY" & plants$short=="CP"] <- opt[3]
  plants$exposed[plants$abbrv=="PVY" & plants$short=="HC-Pro"] <- opt[4]
  plants$exposed[plants$abbrv=="PVY" & plants$short=="NIa-VPg"] <- opt[5]
  plants$exposed[mdat$abbrv=="TMV" & plants$short=="CP"] <- opt[6]
  plants$exposed[mdat$abbrv=="TMV" & plants$short=="RP"] <- opt[7]
  # CP is also the suppressor in ASPV!
  plants$exposed[mdat$abbrv=="ASPV" & plants$short=="CP"] <- opt[8]
  
  prm <- adonis2(wmat ~ log(ncod) + log(treelen) + exposed * enveloped, 
                 data=plants, permutations=9999, by='terms', na.action=na.omit)
  print(prm)
}



####################
##    FIGURE 5    ##
####################


# store mode-specific results
res1 <- adonis2(wmat ~ log(ncod)+log(treelen)+vector, data=mdat, by='terms', permutations=9999)
idx <- mdat$chronic & !mdat$plant
res2 <- adonis2(wmat ~ log(ncod)+log(treelen)+idx, data=mdat, by='terms', permutations=9999)
res3 <- adonis2(wmat ~ log(ncod)+log(treelen)+fecaloral, data=mdat, by='terms', permutations=9999)
res4 <- adonis2(wmat ~ log(ncod)+log(treelen)+respiratory, data=mdat, by='terms', permutations=9999)

pdf(file.path(out.dir, "Figure5.pdf"), width=5, height=5)
par(mar=rep(0.5, 4), mfrow=c(2,2))
plot(m2, type='n', xaxt='n', yaxt='n', bty='n', xlab=NA, ylab=NA)
title(main="Vector-borne", font.main=1, cex.main=1.1, adj=0, line=-0.5)
points(m2, col='grey', cex=ifelse(mdat$vector, 0.01, 0.5), pch=19)
ex <- order(mdat$exposed)
text(m2[ex,], label=labels[ex], cex=ifelse(mdat$vector[ex], 0.6, 0.01), 
     font=ifelse(mdat$exposed[ex], 2, 3), 
     col=ifelse(mdat$exposed[ex], 'darkgreen', rgb(0.05,0.7,0.)))
text(0.99*par('usr')[1], 0.925*par('usr')[3], 
     bquote(R^2 == .(round(100*res1$R2[3], 1))*"%, " ~ P == .(res1$`Pr(>F)`[3])), 
     adj=0, cex=0.8)
abline(v=par('usr')[2], lwd=0.5, xpd=NA)
abline(h=par('usr')[3], lwd=0.5, xpd=NA)
plot(m2, type='n', xaxt='n', yaxt='n', bty='n', xlab=NA, ylab=NA)
title(main="Respiratory", font.main=1, cex.main=1.1, adj=0, line=-0.5)
points(m2, col='grey', cex=ifelse(mdat$respiratory, 0.01, 0.5), pch=19)
text(m2[ex, ], label=labels[ex], xpd=NA,
     cex=ifelse(mdat$respiratory[ex], 0.6, 0.01), 
     font=ifelse(mdat$exposed[ex], 2, 3), 
     col=ifelse(mdat$exposed[ex], 'steelblue4', rgb(0.36,0.54,0.69)))
text(0.99*par('usr')[1], 0.925*par('usr')[3], 
     bquote(R^2 == .(round(100*res4$R2[3], 1))*"%, " ~ P == .(res4$`Pr(>F)`[3])), 
     adj=0, cex=0.8)
plot(m2, type='n', xaxt='n', yaxt='n', bty='n', xlab=NA, ylab=NA)
title(main="Fecal-oral", font.main=1, cex.main=1.1, adj=0, line=-0.5)
points(m2, col='grey', cex=ifelse(mdat$fecaloral, 0.01, 0.5), pch=19)
text(m2[ex,], label=labels[ex], 
     cex=ifelse(mdat$fecaloral[ex], 0.6, 0.01), 
     font=ifelse(mdat$exposed[ex], 2, 3), 
     col=ifelse(mdat$exposed[ex], 'purple4', rgb(0.53, 0.3, 0.74)))
text(0.99*par('usr')[1], 0.925*par('usr')[3], 
     bquote(R^2 == .(round(100*res3$R2[3], 1))*"%, " ~ P == .(res3$`Pr(>F)`[3])), 
     adj=0, cex=0.8)
plot(m2, type='n', xaxt='n', yaxt='n', bty='n', xlab=NA, ylab=NA)
title(main="Chronic STBBI", font.main=1, cex.main=1.1, adj=0, line=-0.5)
points(m2, col='grey', cex=ifelse(idx, 0.01, 0.5), pch=19)
text(m2[ex,], label=labels[ex], cex=ifelse(idx[ex], 0.6, 0.01), 
     font=ifelse(mdat$exposed[ex], 2, 3), 
     col=ifelse(mdat$exposed[ex], 'firebrick4', rgb(0.7, 0.25, 0.25)))
idx <- mdat$chronic & !mdat$plant
text(0.99*par('usr')[1], 0.925*par('usr')[3], 
     bquote(R^2 == .(round(100*res2$R2[3], 1))*"%, " ~ P == .(res2$`Pr(>F)`[3])), 
     adj=0, cex=0.8)
dev.off()


# Supplementary Figure S11
pdf(file.path(out.dir, "mds-contact.pdf"), width=5, height=2.5)
tmodes <- c('vector', 'respiratory', 'fecal-oral', 'STBBI')
pal <- c('forestgreen', 'steelblue', 'purple', 'firebrick')
layout(matrix(c(1,2,5,5,3,4,5,5), nrow=2, byrow=T))
par(mar=c(0,0,0,0))
for (i in 1:4) {
  idx <- (mdat$trans.mode == tmodes[i])
  o <- order(idx, decreasing=F)
  plot(m2[o,], bty='n', xaxt='n', yaxt='n',
       pch=ifelse(idx[o], ifelse(mdat$exposed[o], 19, 1), 19),
       col=ifelse(idx[o], pal[i], 'grey'),
       cex=ifelse(idx[o], 1, 0.3))
  title(main=tmodes[i], adj=0, font.main=1, line=-1, cex.main=0.9)
}
plot(m2, type='n', xaxt='n', yaxt='n', bty='n', xlab=NA, ylab=NA)
abline(v=par('usr')[1], lwd=0.5)
title(main=" Contact, plant", adj=0, font.main=2, line=-1)
points(m2, col='grey', cex=ifelse(mdat$contact, 0.01, 0.5), pch=19)
text(m2, label=labels, cex=ifelse(mdat$contact, 0.8, 0.01), 
     font=2, col='salmon2')
dev.off()

