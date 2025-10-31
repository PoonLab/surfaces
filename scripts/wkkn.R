
# kNN (k = 5), with 10 replicates R = 10

suppressPackageStartupMessages({
  library(caret)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(vegan)
  library(purrr)
})

# Paths 
R <- 10
kval <- 5
in_file <- "/home/hugocastelan/Documents/projects/surfaces_data/L50.RData"          
out_dir   <- "/home/hugocastelan/Documents/projects/ML_kNN_L50"   
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Load  matrix 
load(in_file)  # wmat (NxN), mdat (data.frame)
N      <- nrow(wmat)
n_prot <- N / R

# Ensure symmetry 
ix <- lower.tri(wmat)
wmat[ix] <- t(wmat)[ix]

# Labels 
labels_surfenv <- ifelse(
  mdat$exposed & !is.na(mdat$enveloped),
  ifelse(mdat$enveloped, "enveloped", "non-enveloped"),
  NA
)
labels_surfenv <- factor(labels_surfenv, levels = c("enveloped","non-enveloped"))
cols_env <- c("enveloped" = "#b30000", "non-enveloped" = "#ff8c00")

#Helpers
simple_metrics <- function(conf) {
  acc  <- sum(diag(conf)) / sum(conf)
  prec <- mean(diag(conf) / colSums(conf), na.rm = TRUE)
  rec  <- mean(diag(conf) / rowSums(conf), na.rm = TRUE)
  f1   <- mean(2 * prec * rec / (prec + rec), na.rm = TRUE)
  data.frame(
    Accuracy  = round(100 * acc, 2),
    Precision = round(100 * prec, 2),
    Recall    = round(100 * rec, 2),
    F1        = round(100 * f1, 2)
  )
}

knn_predict_from_D <- function(D_full, train_idx, test_idx, y_train, k = 5) {
  eps <- 1e-9
  preds <- character(length(test_idx))
  for (ii in seq_along(test_idx)) {
    t <- test_idx[ii]
    d <- D_full[t, train_idx]
    ord <- order(d)
    kk <- min(k, length(train_idx))
    nn <- train_idx[ord][seq_len(kk)]
    ww <- 1 / (d[ord][seq_len(kk)] + eps)
    labs <- y_train[match(nn, train_idx)]
    df <- aggregate(ww, by = list(labs), FUN = sum)
    preds[ii] <- as.character(df$Group.1[which.max(df$x)])
  }
  preds
}




# Make a csv 

results_knn <- list()
conf_all <- NULL

for (r in 1:R) {
  idx_rep <- seq(r, N, by = R)
  D_rep   <- wmat[idx_rep, idx_rep, drop = FALSE]
 
  labs_full <- labels_surfenv
  valid     <- which(!is.na(labs_full))
  labs      <- labs_full[valid]
  D_sub     <- D_rep[valid, valid, drop = FALSE]
 
  if (length(unique(labs)) < 2) next
 
  D_sub <- D_sub / max(D_sub)
  diag(D_sub) <- 0
 
  set.seed(123 + r)
  folds <- createFolds(factor(labs), k = 5)
 
  lvls <- levels(factor(labs))
  conf_accum <- matrix(0, nrow = length(lvls), ncol = length(lvls),
                       dimnames = list(lvls, lvls))
  accs <- c(); precs <- c(); recs <- c(); f1s <- c()
 
  for (f in seq_along(folds)) {
    test_idx  <- folds[[f]]
    train_idx <- setdiff(seq_along(labs), test_idx)
    y_train   <- factor(labs[train_idx], levels = lvls)
    y_test    <- factor(labs[test_idx],  levels = lvls)
   
    pred <- knn_predict_from_D(D_sub, train_idx, test_idx, y_train, k = kval)
    conf <- table(Reference = y_test, Prediction = factor(pred, levels = lvls))
   
    mets <- simple_metrics(conf)
    conf_accum <- conf_accum + conf
    accs <- c(accs, mets$Accuracy)
    precs <- c(precs, mets$Precision)
    recs <- c(recs, mets$Recall)
    f1s  <- c(f1s,  mets$F1)
  }
 
  conf_avg <- round(100 * prop.table(conf_accum, margin = 1), 1)
  conf_df <- as.data.frame(as.table(conf_avg))
  colnames(conf_df)[1:3] <- c("Reference", "Prediction", "Freq")
  conf_df$Replicate <- r
  conf_df$k <- kval
  conf_all <- bind_rows(conf_all, conf_df)
 
  results_knn[[length(results_knn) + 1]] <- data.frame(
    Replicate = r, k = kval, Method = "kNN",
    Accuracy = mean(accs), Precision = mean(precs),
    Recall = mean(recs), F1 = mean(f1s)
  )
}

knn_df <- bind_rows(results_knn)
write.csv(knn_df, file.path(out_dir, "kNN_k5_results.csv"), row.names = FALSE)
write.csv(conf_all, file.path(out_dir, "Confusion_k5.csv"),   row.names = FALSE)

valid_reps <- sort(unique(conf_all$Replicate))


#Plot a MSD per replicate 
mds_all <- map_dfr(valid_reps, function(r) {
  idx_rep <- seq(r, N, by = R)
  D_rep   <- wmat[idx_rep, idx_rep, drop = FALSE]
  labs_full <- labels_surfenv
  valid     <- which(!is.na(labs_full))
  labs      <- labs_full[valid]
  D_sub     <- D_rep[valid, valid, drop = FALSE]
  if (length(unique(labs)) < 2) return(NULL)
  D_sub <- D_sub / max(D_sub); diag(D_sub) <- 0
  set.seed(1000 + r)
  mds_fit <- metaMDS(as.dist(D_sub), k = 2, trymax = 50, autotransform = FALSE)
  mds_df  <- as.data.frame(mds_fit$points)
  colnames(mds_df) <- c("MDS1","MDS2")
  mds_df$Label     <- factor(labs, levels = c("enveloped","non-enveloped"))
  mds_df$Replicate <- factor(r, levels = valid_reps)
  mds_df
})

p_mds_row <- ggplot(mds_all, aes(MDS1, MDS2, color = Label)) +
  geom_point(size = 2.8, alpha = 0.95) +
  scale_color_manual(values = cols_env, name = NULL) +
  coord_equal() +
  theme_minimal(base_size = 12) +
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    panel.grid.minor = element_blank()
  ) +
  labs(title = NULL, x = "MDS1", y = "MDS2") +
  facet_wrap(~ Replicate, nrow = 1)   

# Confusion matrix per replicate 

p_conf_row <- conf_all %>%
  mutate(Replicate = factor(Replicate, levels = valid_reps)) %>%
  ggplot(aes(x = Prediction, y = Reference, fill = Freq)) +
  geom_tile(color = "white", linewidth = 0.35) +
  geom_text(aes(label = sprintf("%.1f", Freq)), size = 3.0) +
  scale_fill_gradient(low = "#f7fbff", high = "#08306b") +
  coord_fixed() +
  theme_minimal(base_size = 12) +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    strip.text = element_text(size = 10, face = "bold"),
    panel.grid = element_blank()
  ) +
  labs(title = "kNN (k = 5) – Confusion Matrices per Replicate",
       subtitle = "Rows normalized (%)",
       fill = "Freq") +
  facet_wrap(~ Replicate, nrow = 1)   


# Final plot 

final_plot <- p_mds_row / p_conf_row + plot_layout(heights = c(1, 1.1))


ggsave(file.path(out_dir, "MDS_and_confusion_matrix.png"),
       plot = final_plot, width = 36, height = 14, dpi = 300, bg = "white")
ggsave(file.path(out_dir, "MDS_and_confusion_matrix.pdf"),
       plot = final_plot, width = 36, height = 14, device = cairo_pdf, bg = "white")