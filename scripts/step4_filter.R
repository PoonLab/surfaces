# The purpose of this script is to visualize the decay of tree length
# with the progressive pruning of the shortest tips, and to 
# develop and apply a function for selecting a threshold for pruning.

# Load contents of all CSV files into a single data frame
# produced by running prunetree.py with no additional arguments.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript step4_filter.R \"[glob]\" (PDF output)")
  exit(1)
}
glob <- args[1]
outpath <- ifelse(length(args) > 1, args[2], NA)

files <- Sys.glob(glob)

# Loop through each file to read and process the data
df <- data.frame(protein=character(), ntips=numeric(), tree.len=numeric(), stringsAsFactors=FALSE)
for (fn in files) {
  prune <- read.csv(fn, header=FALSE)
  parts <- strsplit(basename(fn), "_")[[1]]
  prot <- paste(parts[2:(length(parts)-1)], collapse = "_")
  temp <- data.frame(protein=prot, ntips=prune$V1, tree.len=prune$V2, stringsAsFactors=FALSE)
  df <- rbind(df, temp)
}
prots <- unique(df$protein)

#' Locate point where decay slope exceeds some threshold
#' 
#' @param x:  numeric, vector of monotonically decreasing values
#' @param w:  integer, window width (2w+1) to measure slope
#' @param threshold:  numeric, threshold as a fraction of constant decay
#'                    i.e., total tree length / number of tips
#' @return integer, index for point of significant decline
locate.decline <- function(x, w=5, threshold=0.1) {
  # validate input vector
  stopifnot(is.numeric(x))
  stopifnot(length(x) > 2*w+1)  # vector is sufficient length
  stopifnot(all(diff(x) <= 0))  # is monotonically decreasing
  stopifnot(!all(x==x[1]))  # not a flat line

  max.x <- x[1]  # total tree length
  threshold <- threshold * max.x / length(x)  # adjust to constant decay
  
  i <- which(x < max.x)[1]  # fast-forward to initial decline
  while(i < length(x)) {
    slope <- (x[max(1, i-w)] - x[min(i+w, length(x))]) / (2*w+1)
    if (slope > threshold) {
      break
    }
    i <- i+1
  }
  i
}


if (!is.na(outpath)) {
  pdf(outpath, width=5, height=5)
  
  # Visualize tree lengths as function of number of tips
  require(ggfree, quietly = TRUE)
  pal <- ggfree::gg.rainbow(n=length(prots))
  
  # prepare plot region
  par(mar=c(5,5,1,1))
  plot(NA, xlim=range(df$ntips), ylim=range(df$tree.len),
       xlab="Number of tips", ylab="Tree length", bty='n')
}

for (i in 1:length(prots)) {
  prot <- prots[i]
  x <- df$ntips[df$protein==prot]
  y <- df$tree.len[df$protein==prot]
  idx <- locate.decline(y, threshold=0.5)
  print(paste("Protein:", prot, "=", x[idx]))

  if (!is.na(outpath)) {
    lines(x, y, col=pal[i], lwd=2)
    points(x[idx], y[idx], pch=19, col=pal[i])
    text(x[1], y[1], label=prot, col=pal[i], adj=0, cex=0.7, xpd=NA)
    text(0.01*diff(range(df$ntips))+x[idx], y[idx], label=x[idx], cex=0.5, adj=0)
  }
}

if (!is.na(outpath)) {
  dev.off()
}

