set.seed(5)

########################################
# Import k-mer distance
########################################
setwd('/Users/Laura/Projects/surfaces')
km <- read.csv('measles_feb27_kmerdist.csv', header=F, row.names=1)
n.prots <- 8 # number of proteins per genome
km <- as.matrix(km)

########################################
# Hierarchical Clustering
########################################
# t-stochastic neighbour embedding results in more consistent cluster
kd <- 1-km
hc <- hclust(as.dist(kd))
plot(hc, labels=F)
clus <- cutree(hc, h=0.4)
table(clus)

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
                   'clusters'= as.integer(clus))

write.csv(info, 'feb28measles_protein-clusters-info.csv')


