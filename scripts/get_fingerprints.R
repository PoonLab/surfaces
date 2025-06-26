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
  step <- gsub("step([45]).+", "\\1", tokens[length(tokens)])
  protein <- paste(tokens[2:(length(tokens)-1)], collapse=" ")
  
  # posterior distribution
  json <- read_json(f, simplifyVector = T)
  grid <- matrix(json$grid[,3], nrow=20, byrow=T)
  values <- json$grid[1:20, 2]
  
  # empirical distribution
  #mle <- json$MLE$content[["0"]]
  #emp.grid <- fprint(alpha=mle[,1], beta=mle[,2])
  
  list(virus=virus, protein=protein, step=step, values=values, grid=grid) 
       #emp.grid=emp.grid)
}

# default Dirichlet parameter at 0.5
files <- Sys.glob("6_fubar/*.fubar.json")

# load FUBAR grids
grids <- lapply(files, function(f) parse.json(f))
virus <- sapply(grids, function(x) x$virus)
protein <- sapply(grids, function(x) x$protein)
keys <- paste(virus, protein)


mdat <- read.csv("mpp_ap3_metadata.csv", na.strings="")

# match to metadata
mdat$keys <- paste(mdat$virus, mdat$protein)
idx <- which(!is.element(keys, mdat$keys))
length(keys[idx])==0
idx <- which(!is.element(mdat$keys, keys))
idx2 <- which(mdat$include=="N")
all(is.element(idx, idx2))

# take out grid entries that have a 6_fubar/*.json file but should have been
# excluded
remove <- idx2[which(!is.element(idx2, idx))]
mdat[remove, 1:4]
idx <- match(mdat$keys[remove], keys)
grids <- grids[-idx]


write_json(grids, path="fingerprints.json", pretty=TRUE)

# concentration (Dirichlet) parameter set to 1.0
files.2 <- Sys.glob("6_fubar/concentration_1/*.fubar.json")
grids.2 <- lapply(files.2, function(f) parse.json(f))

# visualize results
for (start in seq(1, length(grids), 25)) {
  end <- min(start + 24, length(grids))
  png(paste("fpt", start, "-", end, "_c0.5.png", sep=''), 
      width=6*300, height=6*300, res=300)
  par(mfrow=c(5, 5), mar=c(0,0,0,0))
  for (i in start:end) {
    image(grids[[i]]$grid, xaxt='n', yaxt='n')
    abline(a=0, b=1, col=rgb(0,0,0,0.2))
    text(x=0.5, y=0.95, cex=0.8,
         label=paste(grids[[i]]$virus, grids[[i]]$protein, sep=" "))
  }
  dev.off()  
}


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


