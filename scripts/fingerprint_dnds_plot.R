args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Usage: Rscript \"[glob to fubar CSVs]\" (optional PDF image)")
}
glob <- args[1]
outpath <- ifelse(length(args) > 1, args[2], NA)

fingerprint <- function(
    dnds, 
    breaks = c(0, 0.1, 0.25, 0.5, 1, 1.6, 2.5, 5, 10, 50), 
    cex.lab = 1,
    cex.axis = 0.8, ...
) {
  x <- findInterval(dnds$alpha, vec=breaks)
  y <- findInterval(dnds$beta, vec=breaks)
  tab <- table(x, y)
  n <- length(breaks)
  
  plot(NA, xlim=c(0, n-1), ylim=c(0, n-1), type='n', bty='n',
       xaxt='n', yaxt='n', xlab="Synonymous rates", 
       ylab="Nonsynonymous rates", cex.lab=cex.lab, ...)
  axis(side=1, at=0:(n-1), labels=c(breaks), cex.axis=cex.axis)
  axis(side=2, at=0:(n-1), labels=c(breaks), cex.axis=cex.axis, las=1)
  for (i in 1:dim(tab)[1]) {
    for (j in 1:dim(tab)[2]) {
      z <- (1 - tab[i,j]/max(tab))^3
      if (z<1) {
        z <- min(0.95, z)
      }
      rect(xl=i-1, xr=i, yb=j-1, yt=j, 
           col=rgb(z,z,z), border=NA)
    }
  }
  abline(a=0, b=1, lty=2)
}

# display fingerprints
#setwd("/home/hugocastelan/Documents/projects/surfaces_data/zika1/step5/")
#files <- Sys.glob("/home/hugocastelan/Documents/projects/surfaces_data/zika1/step5/*.fubar.csv")
files <- Sys.glob(glob)
if (!is.na(outpath)) {
  pdf(file=outpath, width=10, height=12)
}

par(mfrow=c(5,4))
for (i in 1:length(files)) {
  filename <- basename(files[i])
  prefix <- strsplit(filename, "\\.")[[1]][1]
  parts <- strsplit(prefix, "_")[[1]]
  prot <- paste(parts[2:(length(parts)-1)], collapse = "_")
  step_no <- as.character(parts[length(parts)])
  method <- ifelse(step_no == "step4", "large", "small")
  dnds <- read.csv(files[i])  
  fingerprint(dnds, main=paste(prot, method))
}

if (!is.na(outpath)) {
  dev.off()
}