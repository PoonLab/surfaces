# The purpose of this script is to visualize the decay of tree length
# with the progressive pruning of the shortest tips, and to 
# develop and apply a function for selecting a threshold for pruning.

# Modify this path to select a specific group of CSV files
# produced by running prunetree.py with no additional arguments.
files <- Sys.glob("~/git/surfaces/data/IAV/*step4.pruning.csv")

# Load contents of all CSV files into a single data frame
prune <- read.csv(files[1], header=F)
df <- data.frame(protein='HA', ntips=prune$V1, tree.len=prune$V2)
for (i in 2:length(files)) {
  prune <- read.csv(files[i], header=F)
  # extract protein name from canonical filename
  prot <- strsplit(basename(files[i]), "_")[[1]][2]
  temp <- data.frame(protein=prot, ntips=prune$V1, tree.len=prune$V2)
  df <- rbind(df, temp)
}


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


# Visualize tree lengths as function of number of tips
require(ggfree)
prots <- unique(df$protein)
pal <- ggfree::gg.rainbow(n=8)

# prepare plot region
par(mar=c(5,5,1,1))
plot(NA, xlim=range(df$ntips), ylim=range(df$tree.len),
     xlab="Number of tips", ylab="Tree length", bty='n')

for (i in 1:length(prots)) {
  prot <- prots[i]
  x <- df$ntips[df$protein==prot]
  y <- df$tree.len[df$protein==prot]
  lines(x, y, col=pal[i], lwd=2)
  
  idx <- locate.decline(y, threshold=1)
  points(x[idx], y[idx], pch=19, col=pal[i])
}



