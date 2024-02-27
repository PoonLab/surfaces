set.seed(5)

########################################
# Import k-mer distance
########################################
setwd('/home/laura/Projects/ProtClust/scripts')
km <- read.csv('measles_distance_filtered.csv', header=F, row.names=1)
n.prots <- 8 # number of proteins per genome
km <- 1-as.matrix(km)

########################################
# t-SNE and Hierarchical Clustering
########################################
# t-stochastic neighbour embedding results in more consistent cluster
# sizes
require(Rtsne)
res <- Rtsne(km, is_distance=T, verbose=T, dims=2)
hc2 <- hclust(dist(res$Y), method='ward.D2')
clusters <- cutree(hc2, n.prots)
#diagnostics.hclust(hc2, x=seq(30, 100, length.out=20))

########################################
# Function to generate color palettes
########################################
gg2.cols <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

########################################
# Plot clustering results
########################################
pal <- gg2.cols(n=max(clusters))
par(mfrow=c(1,1))
plot(res$Y, type='n')
text(res$Y, label=clusters, col=pal[clusters], cex=0.8)

########################################
# Save clustering results as .csv file
########################################
get.part <- function(labels, lab_part=FALSE, sep='-') {
  sapply(labels, function(x) {
    tokens <- strsplit(x, sep)[[1]]  # split the label by dashes
    # tokens[length(tokens)]
    if(lab_part) {tokens[lab_part]}
    else{tokens[length(tokens)]}
    # tokens[lab_part]# return only the last piece
  })
}

# genome labels
labels <- as.character(row.names(km))

# info required to run plot-genome.py
acc <- as.character(get.part(labels, 1))
gene.name <- as.character(get.part(row.names(km), 3))
loc <- as.character(get.part(row.names(km)))

info <- data.frame('name'=row.names(km),
                   'accession'=acc,
                   'coords' = loc,
                   'gene.name' = gene.name,
                   'clusters'= clusters)


write.csv(info, 'measles_protein-clusters-info.csv')


