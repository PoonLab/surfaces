setwd("~/git/surfaces/data")

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
  #step <- gsub("step([45]).+", "\\1", tokens[length(tokens)])
  protein <- paste(tokens[2:(length(tokens)-1)], collapse=" ")
  
  # posterior distribution
  json <- read_json(f, simplifyVector = T)
  grid <- matrix(json$grid[,3], nrow=20, byrow=T)
  values <- json$grid[1:20, 2]
  
  # empirical distribution
  #mle <- json$MLE$content[["0"]]
  #emp.grid <- fprint(alpha=mle[,1], beta=mle[,2])
  
  list(virus=virus, protein=protein, #step=step, 
       values=values, grid=grid) 
       #emp.grid=emp.grid)
}

# default Dirichlet parameter at 0.5
#files <- Sys.glob("6_fubar/*.fubar.json")
#files <- Sys.glob("8_sample/*.fubar.json")  # sample 50 codons
files <- Sys.glob("9_longer/*.fubar.json")  # sample 100 codons

# load FUBAR grids - this takes a minute
grids <- lapply(files, function(f) parse.json(f))

#save(grids, file="step9_grids.RData")  # can use this later in earthmover.R
load("step9_grids.RData")

##################################
# this later code is for visualization of individual fingerprints
grids.0 <- grids[grepl("_0.fubar.json", files)]
grids.1 <- grids[grepl("_1.fubar.json", files)]
grids.2 <- grids[grepl("_2.fubar.json", files)]

virus <- sapply(grids.0, function(x) x$virus)
protein <- sapply(grids.0, function(x) gsub(" step9", "", x$protein))
keys <- paste(virus, protein)


#mdat <- read.csv("mpp_ap3_metadata.csv", na.strings="")
mdat <- read.csv("metadata.csv", na.strings="")

## match to metadata
mdat$keys <- paste(mdat$virus, mdat$protein)
#idx <- which(!is.element(keys, mdat$keys))
#length(keys[idx])==0

idx <- match(keys, mdat$keys)
mdat <- mdat[idx, ]
mdat$label <- paste(mdat$abbrv, mdat$short)
#idx <- which(!is.element(mdat$keys, keys))
#idx2 <- which(mdat$include=="N")
#all(is.element(idx, idx2))


# visualize results
pal <- colorRampPalette(c('white', 'steelblue3', 'grey10'))(24)
for (start in seq(1, length(grids.0), 25)) {
  end <- min(start + 24, length(grids.0))
  #png(paste("fpt.cod50.n0.", start, "-", end, ".png", sep=''), 
  #    width=6*300, height=6*300, res=300)
  pdf(paste("fpt.step9.n0.", start, "-", end, ".pdf", sep=''))
  par(mfrow=c(5, 5), mar=c(0,0,0,0))
  for (i in start:end) {
    image(grids.0[[i]]$grid, xaxt='n', yaxt='n', 
          col=pal)  #hcl.colors(12, 'Purples', rev=T))
    abline(a=0, b=1, col=rgb(0,0,0,0.2))
    text(x=0.5, y=0.95, cex=0.8,
         label=paste(grids.0[[i]]$virus, grids.0[[i]]$protein, sep=" "))
  }
  dev.off()  
}

# do a random sample
set.seed(127)
idx <- sample(1:length(grids.0), 25)
png("~/papers/surfaces/img/fingerprint-sample.pdf", width=6*600, height=6*600, res=600)
par(mar=c(0,0,0,0), mfrow=c(5,5))
for (i in idx) {
  image(grids.0[[i]]$grid, xaxt='n', yaxt='n', 
        col=pal)  #hcl.colors(12, 'Purples', rev=T))
  abline(a=0, b=1, col=rgb(0,0,0,0.2))
  text(x=0.5, y=0.95, cex=1, label=mdat$label[i])
}
dev.off()


# compare replicates
png("set1.png", width=6*300, height=6*300, res=300)
par(mfrow=c(5,5), mar=c(0,0,0,0))
for (i in 1:25) {
  image(grids.1[[i]]$grid, xaxt='n', yaxt='n')
  abline(a=0, b=1, col=rgb(0,0,0,0.2))
  text(x=0.5, y=0.95, cex=0.8, label=keys[i])
}
dev.off()

png("set2.png", width=6*300, height=6*300, res=300)
par(mfrow=c(5,5), mar=c(0,0,0,0))
for (i in 1:25) {
  image(grids.2[[i]]$grid, xaxt='n', yaxt='n')
  abline(a=0, b=1, col=rgb(0,0,0,0.2))
  text(x=0.5, y=0.95, cex=0.8, label=keys[i])
}
dev.off()
