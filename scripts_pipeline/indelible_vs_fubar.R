require(ape)
require(jsonlite)
setwd("/home/laura/Projects/surfaces/scripts_pipeline")

##################################
# Create a random tree
##################################
r.tree <- rtree(100, rooted = TRUE,
                tip.label = NULL)
write.tree(r.tree, file = "random_tree_100.nwk", append = FALSE,
           digits = 10, tree.names = FALSE)

##################################
# Read FUBAR results
##################################
setwd("/home/laura/Projects/surfaces/data_find_tree_length")
files <- Sys.glob("*.FUBAR.json")
names(files) <- files  # so master list has filenames as names
master <- lapply(files, function(f) read_json(f, simplifyVector=TRUE))
names(master)

# Returns MLE estimates for each dataset as list of data frames
# Data frame res[[1]] includes alpha, beta, alpha=beta, LRT, P-value, 
# Total branch length, codon.pos and dnds

res <- lapply(master, function(d) {
  stats <- d$MLE
  site.stats <- as.data.frame(stats$content[[1]][,1:6])
  names(site.stats) <- stats$headers[,1]
  site.stats$codon.pos<-row.names(site.stats)
  site.stats$dnds<- site.stats$beta/site.stats$alpha
  return(site.stats)
})

##################################
# Import indelible run information
##################################

omegas <- c(
  0.0510204, 0.1530612, 0.2551020, 0.3571428, 0.4591836, 0.5612244, 0.6632652, 0.7653060, 0.8673468,
  0.9693876, 1.0714284, 1.1734692, 1.2755100, 1.3775508, 1.4795916, 1.5816324, 1.6836732, 1.7857140,
  1.8877548, 1.9897956, 2.0918364, 2.1938772, 2.2959180, 2.3979588, 2.4999996, 2.6020404, 2.7040812,
  2.8061220, 2.9081628, 3.0102036, 3.1122444, 3.2142852, 3.3163260, 3.4183668, 3.5204076, 3.6224484,
  3.7244892, 3.8265300, 3.9285708, 4.0306116, 4.1326524, 4.2346932, 4.3367340, 4.4387748, 4.5408156,
  4.6428564, 4.7448972, 4.8469380, 4.9489788, 5.0510196)
prop <- c( 
  1.063764e-01, 1.464867e-01, 1.401630e-01, 1.223917e-01, 1.023007e-01, 8.333079e-02, 6.673162e-02, 
  5.279576e-02, 4.139382e-02, 3.222714e-02, 2.495008e-02, 1.922789e-02, 1.476163e-02, 1.129629e-02,
  8.620596e-03, 6.562952e-03, 4.985996e-03, 3.780965e-03, 2.862470e-03, 2.163923e-03, 1.633688e-03,
  1.231906e-03, 9.279287e-04, 6.982665e-04, 5.249686e-04, 3.943506e-04, 2.960034e-04, 2.220245e-04,
  1.664245e-04, 1.246710e-04, 9.333909e-05, 6.984377e-05, 5.223631e-05, 3.904912e-05, 2.917804e-05,
  2.179306e-05, 1.627075e-05, 1.214322e-05, 9.059520e-06, 6.756626e-06, 5.037501e-06, 3.754636e-06,
  2.797655e-06, 2.084012e-06, 1.551998e-06, 1.155507e-06, 8.600995e-07, 6.400653e-07, 4.762155e-07
)

r.files <- Sys.glob("*_RATES.txt")
rates <- lapply(r.files, function(f) read.csv(f, sep="\t"))
rates<- lapply(rates, function(f) {
  # print(length(f$Class))
  f$omegas <- omegas[f$Class+1]
  return(f)
  })
names(rates) <- r.files

#####################################
# Calculate RSME for all tree lengths
#####################################
out <- matrix(ncol = 3, nrow=length(res))

for(i in 1:length(res)){
  name <- names(res)[[i]]
  sel <- res[[i]]
  sim.rates <- rates[[i]]
  
  # --- get tree length from name ---- 
  spl <- strsplit(name, "_")[[1]]
  tree.l <- spl[2]
  
  # --- calculate RSME ---- 
  rsme <- sqrt(mean(sim.rates$omegas - sel$dnds)^2)
  
  # --- store info ---- 
  out[i,1] <-as.numeric(tree.l)
  out[i,2] <-as.numeric(rsme)
  out[i,3] <-as.numeric(spl[3])
} 

out <- as.data.frame(out)
colnames(out) <- c("tree.length", "rmse", "rep")

out <- out[order(out$tree.length),]

#####################################
# Plot RMSE
#####################################

pal <- c("#001219", "#005f73", "#0a9396", 
         "#94d2bd", "#e9d8a6", "#ee9b00",
         "#ca6702", "#bb3e03", "#ae2012",
         "#ae2012")

plot(out$tree.length, out$rmse, 
     col=add.alpha(pal[as.numeric(out$rep)+1], 0.6),
     pch=19)

fit <- loess(rmse ~ tree.length, data = out)
lines(out$tree.length, predict(fit, out),
      col="darkgray", lwd=2, lty=2)

# --- unused---- 
abline(h=0.1, col="gray", lwd=2, lty=2)
abline(v=1, col="gray", lwd=2, lty=2)
  
######################################################
# Unused: Get matching names for rates and fubar files
######################################################
# tree.l <- unlist(
#   regmatches(str,
#              gregexpr("{0,1}[[:digit:]]+\\.{0,1}[[:digit:]]*",str)
#              )
#           )

#   gregexpr("[-]{0,1}[[:digit:]]+\\.{0,1}[[:digit:]]*",str))
# ))