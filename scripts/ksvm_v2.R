# LIBRARIES
library(kernlab)
library(tidyverse)
library(transport)   # Earth Mover's Distance
library(caret)       # Confusion matrix metrics
library(RColorBrewer)
library(patchwork)


# LOAD FILES
load("/Users/hugo/Desktop/step9_grids.RData")   # 'grids'
load("/Users/hugo/Desktop/files.RData")         # 'files'
metadata <- read.csv("/Users/hugo/Desktop/metadata_2.csv", stringsAsFactors = FALSE)


# BUILD FINGERPRINT MATRIX

fp_list <- lapply(grids, function(x) as.vector(x$grid))
lens <- sapply(fp_list, length)
mode_len <- as.numeric(names(sort(table(lens), decreasing = TRUE)[1]))
valid_idx <- which(lens == mode_len)

fp_mat <- do.call(rbind, fp_list[valid_idx])
rownames(fp_mat) <- paste(
  tolower(sapply(grids[valid_idx], `[[`, "virus")),
  tolower(gsub(" step9", "", sapply(grids[valid_idx], `[[`, "protein"))),
  sep = "|"
)

metadata$keys <- paste(tolower(metadata$virus), tolower(metadata$protein), sep = "|")
metadata <- metadata[match(rownames(fp_mat), metadata$keys), ]
dup_idx <- !duplicated(metadata$keys)
fp_mat <- fp_mat[dup_idx, , drop = FALSE]
metadata <- metadata[dup_idx, ]


# NORMALIZE FINGERPRINTS

normalize_row <- function(x) {
  x <- x - min(x)
  if (sum(x) == 0) rep(1 / length(x), length(x)) else x / sum(x)
}
fp_scaled <- t(apply(fp_mat, 1, normalize_row))
rownames(fp_scaled) <- rownames(fp_mat)


# BUILD DISTANCE & KERNEL


#build_matrices <- function(fp_scaled) {
#  n_side <- sqrt(ncol(fp_scaled))
#  coords <- expand.grid(y = 1:n_side, x = 1:n_side)
#  coords <- coords[order(coords$y, coords$x), ]
#  coords <- as.matrix(coords)
  
#  make_wpp <- function(vec) {
#    vec <- as.numeric(vec)
#    if (sum(vec) == 0) vec <- rep(1/length(vec), length(vec))
#    wpp(coords, mass = vec)
#  }
  
#  n <- nrow(fp_scaled)
#  wpps <- lapply(seq_len(n), function(i) make_wpp(fp_scaled[i, ]))
  
#  D <- matrix(0, n, n)
#  for (i in 1:(n - 1)) {
#    for (j in (i + 1):n) {
#      d <- wasserstein(wpps[[i]], wpps[[j]], p = 2, prob = TRUE)
#      D[i, j] <- d; D[j, i] <- d
#    }
#  }
#  Dmax <- max(D)
#  K <- 1 - D / Dmax
#  K[K < 0] <- 0
#  diag(K) <- 1
#  list(D = D, K = K)
#}

#new formula 

build_matrices <- function(fp_scaled, gamma = 1) {
  n_side <- sqrt(ncol(fp_scaled))
  coords <- expand.grid(y = 1:n_side, x = 1:n_side)
  coords <- coords[order(coords$y, coords$x), ]
  coords <- as.matrix(coords)
  
  make_wpp <- function(vec) {
    vec <- as.numeric(vec)
    if (sum(vec) == 0) vec <- rep(1/length(vec), length(vec))
    wpp(coords, mass = vec)
  }
  
  n <- nrow(fp_scaled)
  wpps <- lapply(seq_len(n), function(i) make_wpp(fp_scaled[i, ]))
  
  D <- matrix(0, n, n)
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      d <- wasserstein(wpps[[i]], wpps[[j]], p = 2, prob = TRUE)
      D[i, j] <- d
      D[j, i] <- d
    }
  }
  
  # Change the Kernel by RBF gaussiano simple 
  K <- exp(-gamma * (D^2))
  
  diag(K) <- 1
  return(list(D = D, K = K))
}

# SIMPLE LOO SVM 

run_loo_svm_simple <- function(labels, fp_scaled, title_prefix, Cval = 1) {
  labels <- factor(labels)
  mats <- build_matrices(fp_scaled)
  D <- mats$D; K <- mats$K
  n <- length(labels); preds <- rep(NA, n)
  
  for (i in seq_len(n)) {
    train_idx <- setdiff(seq_len(n), i)
    test_idx  <- i
    svm_fit <- ksvm(
      x = as.kernelMatrix(K[train_idx, train_idx]),
      y = labels[train_idx],
      kernel = "matrix", C = Cval
    )
    preds[i] <- as.character(
      predict(svm_fit, as.kernelMatrix(K[test_idx, train_idx, drop = FALSE]))
    )
  }
  
  conf_mat <- confusionMatrix(factor(preds, levels = levels(labels)), labels)
  
  acc <- round(conf_mat$overall["Accuracy"] * 100, 2)
  prec <- round(conf_mat$byClass["Precision"] * 100, 2)
  rec <- round(conf_mat$byClass["Recall"] * 100, 2)
  f1 <- round(conf_mat$byClass["F1"] * 100, 2)
  
  return(list(D = D, labels = labels, preds = preds,
              conf_mat = conf_mat,
              accuracy = acc, precision = prec, recall = rec, f1 = f1,
              title_prefix = title_prefix))
}


# RUN ANALYSES

results <- list()

# 1. Surface-exposed
valid_idx <- which(!is.na(metadata$exposed))
labels_exposed <- ifelse(metadata$exposed[valid_idx], "positive", "negative")
results$exposed <- run_loo_svm_simple(labels_exposed, fp_scaled[valid_idx, , drop = FALSE],
                                      "Surface-exposed vs Non-exposed")

# 2. Polymerase
valid_idx <- which(!is.na(metadata$polymerase))
labels_poly <- ifelse(metadata$polymerase[valid_idx], "polymerase", "non-polymerase")
results$polymerase <- run_loo_svm_simple(labels_poly, fp_scaled[valid_idx, , drop = FALSE],
                                         "Polymerase vs Non-polymerase")

# 3. Surface vs Polymerase
idx3 <- which(!is.na(metadata$polymerase) & !is.na(metadata$exposed))
labels3 <- ifelse(metadata$polymerase[idx3], "polymerase",
                  ifelse(metadata$exposed[idx3], "surface", NA))
keep3 <- which(!is.na(labels3))
results$surf_vs_poly <- run_loo_svm_simple(labels3[keep3],
                                           fp_scaled[idx3, , drop = FALSE][keep3, , drop = FALSE],
                                           "Surface vs Polymerase")

# 4. Surface-exposed: enveloped vs non-enveloped
idx4 <- which(!is.na(metadata$exposed) & metadata$exposed & !is.na(metadata$enveloped))
labels4 <- ifelse(metadata$enveloped[idx4], "enveloped", "non-enveloped")
results$surf_env <- run_loo_svm_simple(labels4, fp_scaled[idx4, , drop = FALSE],
                                       "Surface-exposed (Enveloped vs Non-enveloped)")

# 5. Polymerase: enveloped vs non-enveloped
idx5 <- which(!is.na(metadata$polymerase) & metadata$polymerase & !is.na(metadata$enveloped))
labels5 <- ifelse(metadata$enveloped[idx5], "enveloped", "non-enveloped")
results$poly_env <- run_loo_svm_simple(labels5, fp_scaled[idx5, , drop = FALSE],
                                       "Polymerase (Enveloped vs Non-enveloped)")


# PLOTS 

pdf("/Users/hugo/Desktop/loo_simple_results.pdf", width = 14, height = 10)

# Plot 1: Surface-exposed vs Non-exposed
mds1 <- cmdscale(results$exposed$D, k = 2)
mds_df1 <- data.frame(X1 = mds1[,1], X2 = mds1[,2],
                      Label = results$exposed$labels,
                      Pred = results$exposed$preds)
p1 <- ggplot(mds_df1, aes(X1, X2, color = Label, shape = Label)) +
  geom_point(size = 4) +
  theme_minimal(base_size = 14) +
  labs(title = paste("MDS:", results$exposed$title_prefix,
                     "| Accuracy:", results$exposed$accuracy, "%"))
cm_df1 <- as.data.frame(results$exposed$conf_mat$table)
colnames(cm_df1) <- c("Real", "Predicted", "Freq")
p1_cm <- ggplot(cm_df1, aes(x = Predicted, y = Real, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Freq), size = 5, fontface = "bold") +
  scale_fill_gradientn(colors = brewer.pal(9, "Blues")) +
  theme_minimal(base_size = 14) +
  labs(title = paste("Confusion Matrix:", results$exposed$title_prefix),
       fill = "Count")
print(p1 + p1_cm)

# Plot 2: Polymerase vs Non-polymerase
mds2 <- cmdscale(results$polymerase$D, k = 2)
mds_df2 <- data.frame(X1 = mds2[,1], X2 = mds2[,2],
                      Label = results$polymerase$labels,
                      Pred = results$polymerase$preds)
p2 <- ggplot(mds_df2, aes(X1, X2, color = Label, shape = Label)) +
  geom_point(size = 4) +
  theme_minimal(base_size = 14) +
  labs(title = paste("MDS:", results$polymerase$title_prefix,
                     "| Accuracy:", results$polymerase$accuracy, "%"))
cm_df2 <- as.data.frame(results$polymerase$conf_mat$table)
colnames(cm_df2) <- c("Real", "Predicted", "Freq")
p2_cm <- ggplot(cm_df2, aes(x = Predicted, y = Real, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Freq), size = 5, fontface = "bold") +
  scale_fill_gradientn(colors = brewer.pal(9, "Blues")) +
  theme_minimal(base_size = 14) +
  labs(title = paste("Confusion Matrix:", results$polymerase$title_prefix),
       fill = "Count")
print(p2 + p2_cm)

# Plot 3: Surface vs Polymerase
mds3 <- cmdscale(results$surf_vs_poly$D, k = 2)
mds_df3 <- data.frame(X1 = mds3[,1], X2 = mds3[,2],
                      Label = results$surf_vs_poly$labels,
                      Pred = results$surf_vs_poly$preds)
p3 <- ggplot(mds_df3, aes(X1, X2, color = Label, shape = Label)) +
  geom_point(size = 4) +
  theme_minimal(base_size = 14) +
  labs(title = paste("MDS:", results$surf_vs_poly$title_prefix,
                     "| Accuracy:", results$surf_vs_poly$accuracy, "%"))
cm_df3 <- as.data.frame(results$surf_vs_poly$conf_mat$table)
colnames(cm_df3) <- c("Real", "Predicted", "Freq")
p3_cm <- ggplot(cm_df3, aes(x = Predicted, y = Real, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Freq), size = 5, fontface = "bold") +
  scale_fill_gradientn(colors = brewer.pal(9, "Blues")) +
  theme_minimal(base_size = 14) +
  labs(title = paste("Confusion Matrix:", results$surf_vs_poly$title_prefix),
       fill = "Count")
print(p3 + p3_cm)

# Plot 4: Surface-exposed (Enveloped vs Non-enveloped)
mds4 <- cmdscale(results$surf_env$D, k = 2)
mds_df4 <- data.frame(X1 = mds4[,1], X2 = mds4[,2],
                      Label = results$surf_env$labels,
                      Pred = results$surf_env$preds)
p4 <- ggplot(mds_df4, aes(X1, X2, color = Label, shape = Label)) +
  geom_point(size = 4) +
  scale_color_manual(values = c("enveloped" = "red",
                                "non-enveloped" = "orange")) +
  theme_minimal(base_size = 14) +
  labs(title = paste("MDS:", results$surf_env$title_prefix,
                     "| Accuracy:", results$surf_env$accuracy, "%"))
cm_df4 <- as.data.frame(results$surf_env$conf_mat$table)
colnames(cm_df4) <- c("Real", "Predicted", "Freq")
p4_cm <- ggplot(cm_df4, aes(x = Predicted, y = Real, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Freq), size = 5, fontface = "bold") +
  scale_fill_gradientn(colors = brewer.pal(9, "Blues")) +
  theme_minimal(base_size = 14) +
  labs(title = paste("Confusion Matrix:", results$surf_env$title_prefix),
       fill = "Count")
print(p4 + p4_cm)

# Plot 5: Polymerase (Enveloped vs Non-enveloped)
mds5 <- cmdscale(results$poly_env$D, k = 2)
mds_df5 <- data.frame(X1 = mds5[,1], X2 = mds5[,2],
                      Label = results$poly_env$labels,
                      Pred = results$poly_env$preds)
p5 <- ggplot(mds_df5, aes(X1, X2, color = Label, shape = Label)) +
  geom_point(size = 4) +
  scale_color_manual(values = c("enveloped" = "red",
                                "non-enveloped" = "orange")) +
  theme_minimal(base_size = 14) +
  labs(title = paste("MDS:", results$poly_env$title_prefix,
                     "| Accuracy:", results$poly_env$accuracy, "%"))
cm_df5 <- as.data.frame(results$poly_env$conf_mat$table)
colnames(cm_df5) <- c("Real", "Predicted", "Freq")
p5_cm <- ggplot(cm_df5, aes(x = Predicted, y = Real, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Freq), size = 5, fontface = "bold") +
  scale_fill_gradientn(colors = brewer.pal(9, "Blues")) +
  theme_minimal(base_size = 14) +
  labs(title = paste("Confusion Matrix:", results$poly_env$title_prefix),
       fill = "Count")
print(p5 + p5_cm)

dev.off()

# SUMMARY CSV
summary_df <- data.frame(
  Analysis = names(results),
  Accuracy = sapply(results, function(x) x$accuracy),
  Precision = sapply(results, function(x) x$precision),
  Recall = sapply(results, function(x) x$recall),
  F1 = sapply(results, function(x) x$f1)
)

write.csv(summary_df, "/Users/hugo/Desktop/loo_simple_summary.csv", row.names = FALSE)
print(summary_df)
