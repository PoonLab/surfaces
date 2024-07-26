# Plot dN (beta) and dS (alpha) from FUBAR analysis
# Analyse FUBAR results for multiple CDSs in a virus
library(plyr)
require(jsonlite)
library(dplyr)

##################################
# Import JSON files 
##################################
setwd("/home/laura/Projects/surfaces/scripts_pipeline/temp_measles/results2")
setwd("/home/laura/Projects/surfaces/evolution_plots/jsons_for_fingerprint_plot")
files <- Sys.glob("*.FUBAR.json")
names(files) <- files  # so master list has filenames as names
master <- lapply(files, function(f) read_json(f, simplifyVector=TRUE))
names(master)

# Returns MLE estimates for each dataset as list of data frames
# Data frame res[[1]] includes alpha, beta, alpha=beta, LRT, P-value, Total branch length, codon.pos and dnds
result <- lapply(master, function(d) {
  stats <- d$MLE
  site.stats <- as.data.frame(stats$content[[1]][,1:6])
  names(site.stats) <- stats$headers[,1]
  site.stats$codon.pos<-row.names(site.stats)
  site.stats$dnds<- site.stats$beta/site.stats$alpha
  return(site.stats)
})

###############################################
# distribution of dN and dS for all proteins
###############################################
par(mfrow=c(3,3)) # 6 rows, 2 columns
# Plot all dn and ds for each protein
for (i in 1:length(result)){
  c1 <- result[[i]]
  prot <- names(result)[[i]]
  header <- gsub(".FUBAR.json", "", prot)
  plot(c1$alpha, main = header, col="dodgerblue", pch=19)
  points(c1$beta, main = header, col="orange", pch=19)
}

###############################################
# Plot a small selection section
###############################################
pdf("/home/laura/Projects/surfaces/evolution_plots/dn_ds_segments.pdf", 
    width=9, height = 3.2)

par(mfrow=c(1,3)) # 6 rows, 2 columns
par(mar=c(5,5,3,2)) # 6 rows, 2 columns
pal <- c("#04447c", "#ff8c00", "#567189", "#A84448")

lims <- list(c(600, 620), c(210, 230), c(420, 440))
n <- c(4, 6, 2)

for (j in 1:length(n)) {
  i <- n[j]
  rdrp <- result[[i]]
  start <- lims[[j]][1]
  end <- lims[[j]][2]
  plot(rdrp$codon.pos, rdrp$alpha, cex=cex.point, pch=21, lwd = 1.3,
       bg=add.alpha(pal[1], 0.7), 
       col=pal[1], main=names(result)[[i]],
       las=1, xlab = "Codon position", ylab = "Rate", type="b", 
       # xlim = c(lims[[i]][1], lims[[i]][2]),
       xlim=c(start, end),
       ylim = c(0,15),
       cex.axis=size.axis, cex.lab=1.5)
  
  points(rdrp$codon.pos, rdrp$beta, cex=cex.point, pch=24, 
         bg=add.alpha(pal[2], 0.7), type = "b", col=pal[2], lwd = 1.3)

}

dev.off()
###############################################
# Plots for Evolution
###############################################
pdf("/home/laura/Projects/surfaces/evolution_plots/rdrp_fingerprint.pdf", 
    width=12, height = 4.2)
par(mfrow=c(1,3)) # 6 rows, 2 columns
par(mar=c(5,5,3,2)) # 6 rows, 2 columns
pal <- c("#04447c", "#ff8c00", "#567189", "#A84448")
i <- 4

# ----- Define values for plots ----- #
cex.point <- 1.7
size.axis <- 1.5
las <- 1

# ----- dN dS distribution ----- #
rdrp <- result[[i]]
# plot(rdrp$alpha, cex=1.5, pch=19, col=add.alpha(pal[1], 0.5), 
#      las=1, xlab = "Codon position", ylab = "Value")
# points(rdrp$beta, cex=1.5, pch=19, col=add.alpha(pal[2], 0.5))
# 

max.y <- max(rdrp$alpha, rdrp$beta)

plot(rdrp$codon.pos, rdrp$alpha, cex=cex.point, pch=21, lwd = 1.3,
     bg=add.alpha(pal[1], 0.7), 
     col=pal[1], main=names(result)[[i]],
     las=1, xlab = "Codon position", ylab = "Rate", type="b", 
     cex.axis=size.axis, cex.lab=1.5)

points(rdrp$codon.pos, rdrp$beta, cex=cex.point, pch=24, 
       bg=add.alpha(pal[2], 0.7), type = "b", col=pal[2], lwd = 1.3)

# ----- Histograms ----- #
br = 100
max.freq <- max(hist(rdrp$alpha, plot=FALSE, breaks=br)$counts,
                      hist(rdrp$beta, plot=FALSE, breaks=br)$counts)

hist(rdrp$alpha, breaks=br, col = add.alpha(pal[1], 0.7),  las = las,
     cex.lab = size.axis, xlab = "Rate", ylab = NA, main="Frequency", 
     ylim = c(0,max.freq), xlim = c(0,10), border = FALSE, cex.axis = size.axis)

hist(rdrp$beta, breaks=br, add=TRUE,col = add.alpha(pal[2], 0.7), border=FALSE)

# ----- Fingerprint ----- #

breaks <- c(0, 0.1, 0.25, 0.5, 1, 1.6, 2.5, 5, 10, 50)
dnds <- rdrp

x <- findInterval(dnds$alpha, vec=breaks)
y <- findInterval(dnds$beta, vec=breaks)
tab <- table(x, y)
n <- length(breaks)

plot(NA, xlim=c(0,n), ylim=c(0,n), type='n', bty='n',
     xaxt='n', yaxt='n', xlab="Synonymous rates", cex.lab=size.axis,
     ylab="Nonsynonymous rates", main=names(result)[[i]])
axis(side=1, at=0:(n-1), labels=c(breaks), cex.axis=size.axis)
axis(side=2, at=0:(n-1), labels=c(breaks), cex.axis=size.axis, las=1)
for (i in 1:dim(tab)[1]) {
  for (j in 1:dim(tab)[2]) {
    z <- (1 - tab[i,j]/max(tab))^3
    if (z<1) {
      z <- min(0.95, z)
    }
    rect(xl=i, xr=i+1, yb=j-1, yt=j, 
         col=rgb(z,z,z), border=NA)
  }
}
abline(a=0, b=1, lty=2)

dev.off()