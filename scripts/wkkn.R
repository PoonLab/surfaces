#require(caret)


#knn <- read.csv("~/papers/surfaces/data/knn.csv", row.names=1)


## Requires wmat and mdatx from earthmover.R
#setwd("~/git/surfaces/data")
load("/Users/hugo/Desktop/L100.RData")  # provides `wmat` and `mdatx`

# k-nearest neighbors
k <- 5
majority <- function(x) {
  vals <- unique(x)
  tab <- tabulate(match(x, vals))
  max.val <- max(tab)
  count <- sum(tab==max.val)
  if (count == 1) {
    return(vals[which.max(tab)])
  } else {
    # tie
    return(NA)
  }
}

pred.stats <- function(true, pred) {
  # assumes binary classification (FALSE, TRUE)
  tab <- table(true, pred)
  tn <- tab[1,1]; fp <- tab[1,2]; fn <- tab[2,1]; tp <- tab[2,2]
  p <- sum(tab[2,]); n <- sum(tab[1,])
  return (c(
    accuracy=(tp+tn)/(p+n),
    precision=tp/(tp+fp),
    f1=2*tp/(2*tp+fp+fn),
    recall=tp/p
  ))
}

################
# surface vs. non-exposed proteins

resolve <- function(votes, weights) {
  totals <- sapply(split(weights, votes), sum)
  as.logical(names(totals)[which.max(totals)])
}

#eps <- 1e-9 #this missing in the original version function

# weighted k-nearest neighbours
wknn <- function(wm, mdat, k=5, eps=1e-9) {
  idx <- which(mdat$exposed)
  pred <- sapply(idx, function(i) {
    # the shortest distance will always be zero (self)
    neighbours <- order(wm[i, idx])[2:(k+1)]
    stopifnot(!is.element(which(mdat$exposed==i), neighbours))
    dists <- wm[i, idx][neighbours]
    weights <- 1/(dists+eps)
    votes <- mdat$enveloped[mdat$exposed][neighbours]
    #majority(votes)
    resolve(votes, weights)
  })
  pred
}

# enveloped vs. non-enveloped surface proteins
# use the first replicate for L=100 samples (step 9)
res <- as.data.frame(t(sapply(1:10, function(j) {
  idx <- seq(j, nrow(wmat), 10)
  mdat <- mdatx[idx, ]
  wm1 <- wmat[idx, idx]  # distances
  
  labels <- mdat$enveloped[mdat$exposed]
  eps <- 1e-9
  pred <- wknn(wm1, mdat)
  
  #table(labels, pred)
  pred.stats(labels, pred)
})))

## show the results 
print(res)

##confusion matrix 
conf.matrices <- lapply(1:10, function(j) {
  idx  <- seq(j, nrow(wmat), 10)
  mdat <- mdatx[idx, ]
  wm1  <- wmat[idx, idx]
  
  labels <- mdat$enveloped[mdat$exposed]
  pred   <- wknn(wm1, mdat)
  
  table(True = labels, Pred = pred)
})

names(conf.matrices) <- paste0("Replicate_", 1:10)

## show the confusion matrix 
conf.matrices



######################################
### Plots
######################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(vegan)
library(purrr)

R <- 10
N <- nrow(wmat)
valid_reps <- 1:10


### Tidy confusion matrices


conf_all <- bind_rows(
  lapply(1:10, function(j) {
    cm <- as.data.frame(conf.matrices[[j]])
    cm$Replicate <- j
    colnames(cm) <- c("Reference", "Prediction", "Freq", "Replicate")
    cm
  })
)


### MDS 
mds_all <- map_dfr(valid_reps, function(r) {
  
  idx_rep <- seq(r, N, by = R)
  D_rep   <- wmat[idx_rep, idx_rep, drop = FALSE]
  
  # FIX: labels must be subset to the replicate
  labs_full <- labels_surfenv[idx_rep] 
  valid     <- which(!is.na(labs_full))
  labs      <- labs_full[valid]
  D_sub     <- D_rep[valid, valid, drop = FALSE]
  
  if (length(unique(labs)) < 2)
    return(NULL)
  
  D_sub <- D_sub / max(D_sub)
  diag(D_sub) <- 0
  
  set.seed(1000 + r)
  mds_fit <- metaMDS(as.dist(D_sub), k = 2, trymax = 50, autotransform = FALSE)
  
  mds_df  <- as.data.frame(mds_fit$points)
  colnames(mds_df) <- c("MDS1","MDS2")
  mds_df$Label     <- factor(labs, levels = c("enveloped","non-enveloped"))
  mds_df$Replicate <- factor(r, levels = valid_reps)
  
  mds_df
})


### MDS – one plot per replicate 

### Colours for the MDS plot


cols_env <- c("enveloped" = "#b30000",
              "non-enveloped" = "#ff8c00")


### MDS Plot 

p_mds_row <- ggplot(mds_all, aes(MDS1, MDS2, color = Label)) +
  geom_point(size = 2.8, alpha = 0.95) +
  scale_color_manual(values = cols_env, name = NULL) +
  coord_equal() +
  theme_minimal(base_size = 12) +
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    panel.grid.minor = element_blank(),
    aspect.ratio = 2  
  ) +
  labs(title = NULL, x = "MDS1", y = "MDS2") +
  facet_wrap(
    ~ Replicate,
    nrow = 1,
    labeller = labeller(Replicate = function(x) paste("Replicate", x))
  )

## Confusion Matrix Plot


p_conf_row <- conf_all %>%
  mutate(
    Reference  = factor(Reference,  levels = c(FALSE, TRUE),
                        labels = c("non-enveloped", "enveloped")),
    Prediction = factor(Prediction, levels = c(FALSE, TRUE),
                        labels = c("non-enveloped", "enveloped")),
    Replicate  = factor(Replicate)
  ) %>%
  ggplot(aes(x = Prediction, y = Reference, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Freq), size = 3.0) +
  scale_fill_gradient(low = "#f7fbff", high = "#08306b") +
  coord_fixed() +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8)
  ) +
  facet_wrap(~ Replicate, nrow = 1,
             labeller = labeller(Replicate = function(x) paste("Replicate", x))) +
  labs(fill = "Freq")

##Text per replice 
res$Replicate <- 1:10

res_lbl <- res %>%
  mutate(
    label = sprintf("Acc = %.3f\nPrec = %.3f\nF1 = %.3f\nRec = %.3f",
                    accuracy, precision, f1, recall),
    Replicate = factor(Replicate)
  ) %>%
  select(Replicate, label)

#panel of text 
p_metrics <- ggplot(res_lbl, aes(x = 1, y = 1, label = label)) +
  geom_text(size = 3.5, lineheight = 1.0) +
  facet_wrap(~ Replicate, nrow = 1) +
  theme_void() +
  theme(
    strip.text = element_text(size = 10, face = "bold")
  )

##Final plot

final_plot <- 
  p_mds_row /
  p_conf_row /
  p_metrics +
  plot_layout(heights = c(1, 1.2, 0.5))

final_plot
