# Cluster homologous sequences from kmer distances

set.seed(5)
########################################
# Modify your filenames
########################################
wd <- '/home/laura/Projects/surfaces'
# input file: .csv file with kmer distances, rownames should be fasta headers
kmer.filename <- 'scripts_pipeline/temp_measles/measles_kmerdist.csv'
info.filename <- 'scripts_pipeline/temp_measles/measles_info2.csv'
########################################
# Import k-mer distance
########################################
setwd(wd)
km <- read.csv(kmer.filename, header=F, row.names=1)
km <- as.matrix(km)
########################################
# Hierarchical Clustering
########################################
# t-stochastic neighbour embedding results in more consistent cluster
kd <- 1-km
hc <- hclust(as.dist(kd))
plot(hc, labels=F)
clus <- cutree(hc, h=0.4)
########################################
# Get protein information
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

# find genome labels
labels <- as.character(row.names(km))

# info required to run plot-genome.py
acc <- as.character(get.part(labels, 1))
gene.name <- as.character(get.part(row.names(km), 3))
loc <- as.character(get.part(row.names(km)))

info <- data.frame('name'=row.names(km),
                   'accession'=acc,
                   'coords' = loc,
                   'gene.name' = gene.name,
                   'clusters'= as.integer(clus),
                   'clus.name'= NA)

# Get cluster name based on most frequent protein annotation
g.n<-table(info$gene.name, info$clusters)
prot.name<-c()
len.clus <-c()
for(i in 1:dim(g.n)[2]){
  m <- which.max(g.n[,i])
  info$clus.name[which(info$clusters==i)] = names(m)
}
########################################
# Save protein information as .csv file
########################################
write.csv(info, info.filename)