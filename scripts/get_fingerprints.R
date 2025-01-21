# generate fingerprints from FUBAR results and export as JSON
# note this is meant to be run on Paphlagon!
setwd("/home/surfaces/")
#setwd("~/git/surfaces/data/")
files <- Sys.glob("6_fubar/*.fubar.csv")

parse.file <- function(f) {
  tokens <- strsplit(basename(f), "_")[[1]]
  virus <- tokens[1]
  step <- gsub("step([45]).+", "\\1", tokens[length(tokens)])
  protein <- paste(tokens[2:(length(tokens)-1)], collapse=" ")
  dnds <- read.csv(f)
  data.frame(virus=virus, protein=protein, step=step,
             pos=dnds$pos, alpha=dnds$alpha, beta=dnds$beta)
}

fubar <- do.call(rbind, lapply(files, parse.file))

manual <- TRUE
nb <- 7
a.breaks <- 10^quantile(log10(info$alpha), prob=seq(0, 1, length.out=nb))
b.breaks <- 10^quantile(log10(info$beta), prob=seq(0, 1, length.out=nb))
if (manual) {
  a.breaks <- 10^c(-3, -0.2, 0, 0.2, 0.4, 0.5, 0.7, 0.85, 2)
  b.breaks <- 10^c(-5, -1.9, -1.5, -1.2, -0.9, -0.5, -0.2, 0.2, 2)
}

fprint <- function(dnds) {
  x <- findInterval(dnds$alpha, vec=a.breaks)
  y <- findInterval(dnds$beta, vec=b.breaks)
  xx <- factor(x, levels=1:(length(a.breaks)-1))
  yy <- factor(y, levels=1:(length(b.breaks)-1))
  tab <- table(xx, yy)
  matrix(tab/sum(tab), ncol=ncol(tab), dimnames=dimnames(tab))
}

by.gene <- lapply(
  split(fubar, fubar$virus), function(df) {
    lapply(split(df, df$protein), function(df2) {
      step4 <- NA
      if (any(df2$step==4)) step4 <- fprint(df2[df2$step==4,])
      step5 <- NA
      if (any(df2$step==5)) step5 <- fprint(df2[df2$step==5,])
      list('step4'=step4, 'step5'=step5)
    })
})

write_json(by.gene, "fingerprints.json", pretty=TRUE)
