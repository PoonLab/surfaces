setwd("~/git/surfaces/data/")

fubar <- read.csv("march21_fubar.csv")
fubar <- fubar[fubar$tree.length >= 0.1, ]

fel <- read.csv("march21_fel.csv")
fel <- fel[fel$tree.length >= 0.1, ]

mdat <- read.csv("iss135/metadata_iss135.csv")
mdat$key <- paste(mdat$virus, mdat$protein)
idx <- match(keys, mdat$key)
mdat <- mdat[idx, ]
astats <- read.csv("iss135/iss135_align_stats.csv")
astats$key <- paste(astats$virus, astats$protein)
idx <- match(mdat$key, astats$key)
mdat$nseq <- astats$nseq[idx]
mdat$ncod <- astats$ncod[idx]
mdat$treelen <- astats$treelen[idx]


pdf("~/papers/surfaces/img/RMSE-tree-length.pdf", width=9, height=4)
par(mfrow=c(1,2), mar=c(5,5,1,5))
plot(jitter(fubar$tree.length), fubar$rmse, bty='n', log='x', pch=19, cex=0.5, 
     col=rgb(0,0,0,0.3), xlab="Tree length (codons)", ylab="RMSE", main='FUBAR', 
     las=2, font.main=1)
rect(xl=1.5, xr=6, yb=0, yt=1.5, col=rgb(0,0,1,0.1), border=NA)
y <- sapply(split(fubar$rmse, fubar$tree.length), mean)
x <- as.numeric(names(y))
lines(names(y), y)
points(names(y), y, pch=21, cex=1.5, col='white', bg='black')

at.y <- seq(par('usr')[3], par('usr')[4], length.out=5)
axis(side=4, at=at.y, labels=seq(0, 1, length.out=5), 
     col='red', col.axis='red', las=2)
mtext("Trees below threshold", side=4, line=3,col='red')
y <- sapply(x, function(tl) sum(3*mdat$treelen < tl) / nrow(mdat))
dy <- diff(par('usr')[3:4])
yy <- y*dy + par('usr')[3]
lines(x, yy, col='red', lty=2)


plot(jitter(fel$tree.length), fel$rmse, bty='n', log='xy', pch=19, cex=0.5, 
     col=rgb(0,0,0,0.3), xlab="Tree length (codons)", ylab="log(RMSE)", 
     main='FEL', font.main=1, las=2)
rect(xl=1.5, xr=6, yb=.1, yt=1e4, col=rgb(0,0,1,0.1), border=NA)
y <- sapply(split(fel$rmse, fel$tree.length), mean)
x <- as.numeric(names(y))
lines(names(y), y)
points(names(y), y, pch=21, cex=1.5, col='white', bg='black')

at.y <- seq(par('usr')[3], par('usr')[4], length.out=5)
axis(side=4, at=10^(at.y), labels=seq(0, 1, length.out=5), 
     col='red', col.axis='red', las=2)
mtext("Trees below threshold", side=4, line=3,col='red')
y <- sapply(x, function(tl) sum(3*mdat$treelen < tl) / nrow(mdat))
dy <- diff(par('usr')[3:4])
yy <- y*dy + par('usr')[3]
lines(x, 10^yy, col='red', lty=2)

dev.off()