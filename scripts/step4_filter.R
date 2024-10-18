#!/usr/bin/env Rscript

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if directory path is provided
if (length(args) < 1) {
  stop("Usage: Rscript step4_filter.R [\"PATH/*_step4.pruning.csv\"] (optional PNG)")
}

# Check if output path is provided (optional)
outpath <- ifelse(length(args) > 1, args[2], NA)

# Load contents of all CSV files into a single data frame
files <- Sys.glob(args[1])

# Check if any CSV files were found
if (length(files) == 0) {
  stop("No CSV files found in the specified directory.")
}

# Load required library
require(ggfree, quietly=TRUE)

# Initialize an empty data frame
df <- data.frame(protein=character(), ntips=numeric(), tree.len=numeric(), stringsAsFactors=FALSE)

# Loop through each file to read and process the data
for (i in 1:length(files)) {
  prune <- read.csv(files[i], header=FALSE)
  filename <- basename(files[i])
  parts <- strsplit(filename, "_")[[1]]
  prot <- paste(parts[2:(length(parts)-1)], collapse = "_")
  temp <- data.frame(protein=prot, ntips=prune$V1, tree.len=prune$V2, stringsAsFactors=FALSE)
  df <- rbind(df, temp)
}

# Function to locate point of significant tree length decline
locate.decline <- function(x, w=5, threshold=0.1) {
  stopifnot(is.numeric(x))
  stopifnot(length(x) > 2*w+1)
  stopifnot(all(diff(x) <= 0))
  stopifnot(!all(x == x[1]))
  
  max.x <- x[1]
  threshold <- threshold * max.x / length(x)
  
  i <- which(x < max.x)[1]
  while (i < length(x)) {
    slope <- (x[max(1, i-w)] - x[min(i+w, length(x))]) / (2*w+1)
    if (slope > threshold) {
      break
    }
    i <- i + 1
  }
  return(i)
}

# Visualize tree lengths as a function of number of tips
prots <- unique(df$protein)
pal <- ggfree::gg.rainbow(n=length(prots))

# Prepare plot region and generate the plot
par(mar=c(3, 3, 2, 2))
plot(NA, xlim=range(df$ntips), ylim=range(df$tree.len), 
     xlab="Number of tips", ylab="Tree length", bty='n')

for (i in 1:length(prots)) {
  prot <- prots[i]
  x <- df$ntips[df$protein == prot]
  y <- df$tree.len[df$protein == prot]
  lines(x, y, col=pal[i], lwd=2)
  
  idx <- locate.decline(y, threshold=0.5)
  points(x[idx], y[idx], pch=19, col=pal[i])
  text(x[1], y[1], label=prot, col=pal[i], adj=0, cex=0.7, xpd=NA)
  text(0.01 * diff(range(df$ntips)) + x[idx], y[idx], label=x[idx], cex=0.5, adj=0)
  print(paste("Protein:", prot, "=", x[idx]))
}

# Save the plot as a PDF if an output path is provided
if (!is.na(outpath)) {
  png(outpath, width=5*150, height=5*150, res=150)
  par(mar=c(5,5,1,1))
  plot(NA, xlim=range(df$ntips), ylim=range(df$tree.len), 
       xlab="Number of tips", ylab="Tree length", bty='n')
  
  for (i in 1:length(prots)) {
    prot <- prots[i]
    x <- df$ntips[df$protein == prot]
    y <- df$tree.len[df$protein == prot]
    lines(x, y, col=pal[i], lwd=2)
    
    idx <- locate.decline(y, threshold=0.5)
    points(x[idx], y[idx], pch=19, col=pal[i])
    text(x[1], y[1], label=prot, col=pal[i], adj=0, cex=0.7, xpd=NA)
    text(0.01 * diff(range(df$ntips)) + x[idx], y[idx], label=x[idx], cex=0.5, adj=0)
  }
  dev.off()
}
