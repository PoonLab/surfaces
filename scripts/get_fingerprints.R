breaks <- (50*(0:19)^5)/19^5  # from Murrell et al. 2016

# generate fingerprints from FUBAR results and export as JSON
fprint <- function(alpha, beta, breaks=breaks) {
  x <- findInterval(alpha, vec=breaks)
  y <- findInterval(beta, vec=breaks)
  xx <- factor(x, levels=1:(length(breaks)-1))
  yy <- factor(y, levels=1:(length(breaks)-1))
  tab <- table(xx, yy)
  matrix(tab/sum(tab), ncol=ncol(tab), dimnames=dimnames(tab))
}


# read FUBAR JSON and concatenate posterior surfaces
require(jsonlite)

parse.json <- function(f) {
  # extract information from filename
  tokens <- strsplit(basename(f), "_")[[1]]
  virus <- tokens[1]
  step <- gsub("step([45]).+", "\\1", tokens[length(tokens)])
  protein <- paste(tokens[2:(length(tokens)-1)], collapse=" ")
  
  # posterior distribution
  json <- read_json(f, simplifyVector = T)
  grid <- matrix(json$grid[,3], nrow=20, byrow=T)
  values <- json$grid[1:20, 2]
  
  # empirical distribution
  mle <- json$MLE$content[["0"]]
  emp.grid <- fprint(alpha=mle[,1], beta=mle[,2])
  
  list(virus=virus, protein=protein, step=step, values=values, grid=grid, 
       emp.grid=emp.grid)
}

# default Dirichlet parameter at 0.5
files <- Sys.glob("6_fubar/*.fubar.json")
grids <- lapply(files, function(f) parse.json(f))

write_json(grids, path="fingerprints.json", pretty=TRUE)

# concentration (Dirichlet) parameter set to 1.0
files.2 <- Sys.glob("6_fubar/concentration_1/*.fubar.json")
grids.2 <- lapply(files.2, function(f) parse.json(f))

# visualize results
start <- 51
end <- start + 24
png(paste("fpt", start, "-", end, "_c0.5.png", sep=''), 
    width=6*300, height=6*300, res=300)
par(mfrow=c(5, 5), mar=c(0,0,0,0))
for (i in start:end) {
  image(grids[[i]]$emp.grid, xaxt='n', yaxt='n')
  abline(a=0, b=1, col=rgb(0,0,0,0.2))
  text(x=0.5, y=0.95, cex=0.8,
       label=paste(grids[[i]]$virus, grids[[i]]$protein, sep=" "))
}
dev.off()

png(paste("fpt", start, "-", end, "_c1.0.png", sep=''), 
    width=6*300, height=6*300, res=300)
par(mfrow=c(5, 5), mar=c(0,0,0,0))
for (i in start:end) {
  image(grids.2[[i]]$grid, xaxt='n', yaxt='n')
  abline(a=0, b=1, col=rgb(0,0,0,0.2))
  text(x=0.5, y=0.95, cex=0.8,
       label=paste(grids.2[[i]]$virus, grids.2[[i]]$protein, sep=" "))
 }
dev.off()


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
write.csv(wmat, file="wdist.csv", quote=F)

wdist <- as.dist(wmat)

mds <- cmdscale(wdist, k=2)
par(mfrow=c(1,1))
plot(mds[,1:2], type='n')
text(mds[,1], mds[,2], labels=names(wpps), cex=0.5)


hc <- hclust(wdist)
pdf("hclust.pdf", width=15, height=5)
par(mar=c(0,0,0,0))
plot(hc, cex=0.5, main=NA)
dev.off()
