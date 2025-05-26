library(jsonlite)
library(kernlab)
library(tidyverse)

setwd("~/git/surfaces/data")

# load json
json_data <- fromJSON("fingerprints.json")
flatten_matrix <- function(mat) {
  as.vector(t(mat))
}
feature_list <- list()
row_names <- c()
protein_names <- c()
virus_names <- c()
for (virus in names(json_data)) {
  for (protein in names(json_data[[virus]])) {
    mat <- json_data[[virus]][[protein]]$step5
    if (!is.null(mat)) {
      vec <- flatten_matrix(mat)
      full_name <- paste0("virus_", virus, "_", protein)
      feature_list[[full_name]] <- vec
      row_names <- c(row_names, full_name)
      protein_names <- c(protein_names, protein)
      virus_names <- c(virus_names, virus)
    }
  }
}

df <- do.call(rbind, feature_list)
df <- as.data.frame(df)
rownames(df) <- row_names
df$protein <- protein_names
df$virus <- virus_names

# load metadata
#annotation <- read.csv("/Users/hugo/Downloads/annotation_with_MPP_AP3.csv")
annotation <- read.csv("annotation_with_MPP_AP3.csv")
df_labeled <- merge(df, annotation, by = "protein")

#Normalization

# Classification labels
labels <- as.factor(df_labeled$Stability)

# Extract numeric columns (fingerprints only)
feature_cols <- sapply(df_labeled, is.numeric)
fingerprint_matrix <- df_labeled[, feature_cols]

# Filter valid rows (no NAs)
valid_rows <- complete.cases(fingerprint_matrix)
fingerprint_clean <- fingerprint_matrix[valid_rows, ]
labels_clean <- labels[valid_rows]

# Normalize
fingerprint_scaled <- scale(fingerprint_clean)

# kernel PCA
sim_matrix <- as.matrix(fingerprint_scaled) %*% t(as.matrix(fingerprint_scaled))
sim_matrix[!is.finite(sim_matrix)] <- 0
km <- as.kernelMatrix(sim_matrix)

# Run KPCA
p <- kpca(km)

# Explained variance (%)
eigenvalues <- eig(p)
var_exp <- eigenvalues / sum(eigenvalues)
pc1_percent <- round(var_exp[1] * 100, 1)
pc2_percent <- round(var_exp[2] * 100, 1)

# Coordinates
coords <- rotated(p)[, 1:2]

# plot KPCA
label_colors <- hcl.colors(2)
plot(coords,
     col = label_colors[labels_clean],
     pch = 19,
     xlab = paste0("PC1 (", pc1_percent, "%)"),
     ylab = paste0("PC2 (", pc2_percent, "%)"),
     main = "KPCA: Classification by Stability")
legend("topright", legend = levels(labels_clean),
       col = label_colors, pch = 19)

#Train and evaluate
fit <- ksvm(as.matrix(fingerprint_scaled), labels_clean, kernel = "rbfdot", C = 1)
pred <- predict(fit)
cat("\nConfusion Matrix:\n")
print(table(Real = labels_clean, Predicted = pred))
cat("\nSupport Vectors:\n")
print(SVindex(fit))



# extract fingerprints from JSON, choosing step4 when step5 is not available
# run earth mover's distance on fingerprints to get distance matrix
# convert this distance matrix into a kernel matrix (similarity = 1-distance)

# train and test
n <- length(labels_clean)
train <- sample(1:n, n/2)  # random subset of 50% of data
mx <- as.matrix(fingerprint_scaled)
fit <- ksvm(mx[train, train], labels_clean[train], kernel="rbfdot", )

