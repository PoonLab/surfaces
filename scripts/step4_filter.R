# The purpose of this script is to visualize the decay of tree length
# with the progressive pruning of the shortest tips, and to 
# develop and apply a function for selecting a threshold for pruning.

# Load contents of all CSV files into a single data frame
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Usage: Rscript \"[glob to pruning CSVs]\" (optional PDF image)")
}

glob <- args[1]
outpath <- ifelse(length(args) > 1, args[2], NA)

files <- Sys.glob(glob)

# Loop through each file to read and process the data
df <- data.frame(protein = character(), ntips = numeric(), tree.len = numeric(), stringsAsFactors = FALSE)

for (fn in files) {
  prune <- read.csv(fn, header = FALSE)
  parts <- strsplit(basename(fn), "_")[[1]]
  prot <- paste(parts[2:(length(parts)-1)], collapse = "_")
  temp <- data.frame(protein = prot, ntips = prune$V1, tree.len = prune$V2, stringsAsFactors = FALSE)
  df <- rbind(df, temp)
}
prots <- unique(df$protein)

# Function to locate point where decay slope exceeds some threshold
locate.decline <- function(x, w = 5, threshold = 0.1) {
  # Validate input vector
  stopifnot(is.numeric(x))
  stopifnot(length(x) > 2 * w + 1)  # Vector is sufficient length
  stopifnot(all(diff(x) <= 0))  # Is monotonically decreasing
  stopifnot(!all(x == x[1]))  # Not a flat line
  
  max.x <- x[1]  # Total tree length
  threshold <- threshold * max.x / length(x)  # Adjust to constant decay
  
  i <- which(x < max.x)[1]  # Fast-forward to initial decline
  while (i < length(x)) {
    slope <- (x[max(1, i - w)] - x[min(i + w, length(x))]) / (2 * w + 1)
    if (slope > threshold) {
      break
    }
    i <- i + 1
  }
  i
}

<<<<<<< HEAD
# Define prots outside of the conditional to avoid 'not found' error
=======
>>>>>>> cc06f3932d17f9797b9f055170558ca04dcd07e1
prots <- unique(df$protein)

if (!is.na(outpath)) {
  pdf(outpath, width = 5, height = 5)
  
  # Visualize tree lengths as function of number of tips
  require(ggfree, quietly = TRUE)
<<<<<<< HEAD
  pal <- ggfree::gg.rainbow(n = length(prots))
=======
  pal <- ggfree::gg.rainbow(n=length(prots))
>>>>>>> cc06f3932d17f9797b9f055170558ca04dcd07e1
  
  # Prepare plot region
  par(mar = c(5, 5, 1, 1))
  plot(NA, xlim = range(df$ntips), ylim = range(df$tree.len),
       xlab = "Number of tips", ylab = "Tree length", bty = 'n')
}

# Loop through proteins to plot data
for (i in seq_along(prots)) {
  prot <- prots[i]
  x <- df$ntips[df$protein == prot]
  y <- df$tree.len[df$protein == prot]
  idx <- locate.decline(y, threshold = 0.5)
  print(paste("Protein:", prot, "=", x[idx]))
  
  if (!is.na(outpath)) {
    lines(x, y, col = pal[i], lwd = 2)
    points(x[idx], y[idx], pch = 19, col = pal[i])
    text(x[1], y[1], label = prot, col = pal[i], adj = 0, cex = 0.7, xpd = NA)
    text(0.01 * diff(range(df$ntips)) + x[idx], y[idx], label = x[idx], cex = 0.5, adj = 0)
  }
}

if (!is.na(outpath)) {
  dev.off()
}
