# Load libraries
library(jsonlite)     
library(kernlab)       # For SVM
library(tidyverse)     
library(transport)     # For computing EMD distances


# Function to flatten fingerprint matrices
flatten_matrix <- function(mat) {
  if (is.null(mat)) return(NULL)
  mat <- as.matrix(mat)
  if (!is.matrix(mat) || any(dim(mat) != c(8, 8))) return(NULL)
  as.vector(t(mat))  
}

# Load and process JSON data
json_data <- fromJSON("/Users/hugo/fingerprints.json")

# Initialize storage structures
feature_list <- list()
row_names <- c()
protein_names <- c()
virus_names <- c()
discarded <- 0

# Process JSON entries
for (virus in names(json_data)) {
  for (protein in names(json_data[[virus]])) {
    entry <- json_data[[virus]][[protein]]
    mat <- entry$step5 %||% entry$step4
    vec <- flatten_matrix(mat)
    
    if (is.null(vec) || any(!is.finite(vec)) || length(vec) != 64) {
      discarded <- discarded + 1
      next
    }
    
    id <- paste0("virus_", virus, "_", protein)
    feature_list[[id]] <- vec
    row_names <- c(row_names, id)
    protein_names <- c(protein_names, paste0(virus, "|", protein))
    virus_names <- c(virus_names, virus)
  }
}

cat("Entries processed:", length(feature_list), "\nDiscarded:", discarded, "\n")

# Build main dataframe
df <- as.data.frame(do.call(rbind, feature_list))
rownames(df) <- row_names
df$protein <- protein_names
df$virus <- virus_names



# Load and process annotations
annotation <- read.csv("/Users/hugo/Downloads/annotation_with_MPP_AP3.csv")
annotation$protein <- paste0(annotation$virus, "|", annotation$protein)

# Filter and consolidate annotations
annotation_filtered <- annotation %>%
  filter(protein %in% df$protein) %>%
  group_by(protein) %>%
  summarize(Stability = names(which.max(table(Stability))), .groups = 'drop')

# Join with main data
df_labeled <- inner_join(df, annotation_filtered, by = "protein")
labels <- factor(df_labeled$Stability)


# Row normalization function
normalize_row <- function(x) {
  x <- x - min(x)
  if (sum(x) == 0) rep(1 / length(x), length(x)) else x / sum(x)
}

# Normalize fingerprint vectors
fingerprint_cols <- which(sapply(df_labeled, is.numeric))[1:64]
fp_clean <- df_labeled[, fingerprint_cols]
valid <- complete.cases(fp_clean)
fp_clean <- fp_clean[valid, ]
labels_clean <- labels[valid]
fp_scaled <- t(apply(fp_clean, 1, normalize_row))

# Prepare coordinates for WPP
coords <- expand.grid(y = 1:8, x = 1:8)
coords <- coords[order(coords$y, coords$x), ]
coords <- as.matrix(coords)

#run earth mover's distance on fingerprints to get distance matrix
make_wpp <- function(vec) {
  vec <- as.numeric(vec)
  vec <- vec / sum(vec)
  wpp(coords, mass = vec)
}

# Create WPP objects
wpps <- lapply(seq_len(nrow(fp_scaled)), function(i) make_wpp(fp_scaled[i, ]))

# Compute EMD distance matrix
n <- length(wpps)
D <- matrix(0, n, n)

for (i in 1:(n - 1)) {
  for (j in (i + 1):n) {
    d <- wasserstein(wpps[[i]], wpps[[j]], p = 2, prob = TRUE)
    D[i, j] <- d
    D[j, i] <- d
  }
}

# MACHINE LEARNING MODELING
#convert this distance matrix into a kernel matrix (similarity = 1-distance)
# Create kernel matrix
K <- 1 - D
# Ensure values are in [0, 1] 
K[K < 0] <- 0


# Train/test split
set.seed(123) #reproducibility data 
train_idx <- sample(seq_len(n), floor(0.7 * n))  # 70% training
test_idx <- setdiff(seq_len(n), train_idx)       # 30% test

# Train SVM model
svm_fit <- ksvm(
  x = as.kernelMatrix(K[train_idx, train_idx]),
  y = labels_clean[train_idx],
  kernel = "matrix",
  C = 1
)

# Evaluate model
pred <- predict(svm_fit, as.kernelMatrix(K[test_idx, train_idx]))
conf_mat <- table(Real = labels_clean[test_idx], Predicted = pred)
print(conf_mat)

# Compute accuracy
accuracy <- sum(diag(conf_mat)) / sum(conf_mat)
cat(sprintf("Accuracy: %.2f%%\n", accuracy * 100))



##VISUALIZATION PLOT
pcoa_result <- pcoa(D)
plot_data <- tibble(
      PC1 = pcoa_result$vectors[, 1],
      PC2 = pcoa_result$vectors[, 2],
      Stability = labels_clean[valid]
)

ggplot(plot_data, aes(x = PC1, y = PC2, color = Stability)) +
     geom_point(size = 3, alpha = 0.8) +
     stat_ellipse(level = 0.95) +
     scale_color_viridis_d() +
     labs(
          title = "PCoA base on EMD",
          x = paste0("PC1 (", round(pcoa_result$values$Relative_eig[1] * 100, 1), "%)"),
          y = paste0("PC2 (", round(pcoa_result$values$Relative_eig[2] * 100, 1), "%)")
       ) +  theme_minimal()
