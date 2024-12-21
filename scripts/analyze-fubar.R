fingerprint <- function(
    dnds, 
    breaks = c(0, 0.1, 0.25, 0.5, 1, 1.6, 2.5, 5, 10, 50)
) {
  x <- findInterval(dnds$alpha, vec=breaks)
  y <- findInterval(dnds$beta, vec=breaks)
  xx <- factor(x, levels=1:10)
  yy <- factor(y, levels=1:10)
  tab <- table(xx, yy)
  tab/sum(tab)  # normalize
}

# from Hugo
cosine <- function(x, y) {
  1 - sum(x * y) / sqrt(sum(x * x) * sum(y * y))  
}

setwd("~/git/surfaces")
fubar <- read.csv("data/collated-fubar.csv")
fubar <- fubar[fubar$step==5, ]

by.gene <- split(fubar, f=list(fubar$protein, fubar$virus), drop=TRUE)
prints <- lapply(by.gene, fingerprint)

# cosine distances
n <- length(prints)
csd <- matrix(0, nrow=n, ncol=n)
for (i in 1:(n-1)) {
  for (j in (i+1):n) {
    csd[i,j] <- csd[j,i] <- cosine(prints[[i]], prints[[j]])
  }
}


# visualize
mds <- cmdscale(csd)
par(mfrow=c(1,1), mar=c(0,0,0,0))
plot(mds, type='n')
text(mds[,1], mds[,2], label=names(by.gene), cex=0.4)

hc <- hclust(as.dist(csd))
plot(hc)

