# -------------------------- Load packages --------------------------
library(openxlsx)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(SingleCellExperiment)
library(scater)

# -------------------------- Files (place in working directory) --------------------------
count_file <- "matrix.xlsx"
group_file <- "class_info_DESeq2.xlsx"
output_dir <- "UMAP_Results"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# -------------------------- Read expression matrix (FPKM) --------------------------
cat("Reading FPKM matrix:", count_file, "\n")
df_raw <- read.xlsx(count_file, rowNames = FALSE, check.names = FALSE)

gene_names <- as.character(df_raw[, 1])

if (any(is.na(gene_names))) {
  na_count <- sum(is.na(gene_names))
  cat("Found", na_count, "rows with missing gene names, removing them\n")
  keep <- !is.na(gene_names)
  df_raw <- df_raw[keep, ]
  gene_names <- gene_names[keep]
}

dup_flag <- duplicated(gene_names)
if (any(dup_flag)) {
  cat("Found", sum(dup_flag), "duplicate gene names, keeping first occurrence\n")
  df_raw <- df_raw[!dup_flag, ]
  gene_names <- gene_names[!dup_flag]
}

rownames(df_raw) <- gene_names
fpkm_data <- as.matrix(df_raw[, -1, drop = FALSE])
fpkm_data[is.na(fpkm_data)] <- 0

n_samples <- ncol(fpkm_data)
cat("Number of samples:", n_samples, "\n")

min_samples <- 27
keep_genes <- rowSums(fpkm_data > 0) >= min_samples
fpkm_data <- fpkm_data[keep_genes, ]
cat("Genes retained after filtering:", nrow(fpkm_data), "\n")

if (nrow(fpkm_data) == 0) stop("No genes retained, lower min_samples or check expression matrix")

sample_names <- colnames(fpkm_data)
sample_names <- trimws(sample_names)
sample_names <- gsub("-", ".", sample_names)
colnames(fpkm_data) <- sample_names

# -------------------------- Read grouping information --------------------------
cat("Reading grouping file:", group_file, "\n")
group_df <- read.xlsx(group_file)
colnames(group_df) <- c("sample", "group")
group_df$sample <- trimws(group_df$sample)
group_df$sample <- gsub("-", ".", group_df$sample)

if (!all(sample_names %in% group_df$sample)) {
  missing <- setdiff(sample_names, group_df$sample)
  cat("Samples missing in grouping file:", paste(missing, collapse = ", "), "\n")
  stop("Sample names in expression matrix and grouping file do not match")
}
group_df <- group_df[match(sample_names, group_df$sample), , drop = FALSE]

cell_meta <- data.frame(
  sample = sample_names,
  group = group_df$group,
  row.names = sample_names,
  stringsAsFactors = FALSE
)
cat("Group summary:\n")
print(table(cell_meta$group))

# -------------------------- Log transformation --------------------------
fpkm_range <- range(fpkm_data, na.rm = TRUE)
cat("FPKM range:", fpkm_range, "\n")
if (max(fpkm_data, na.rm = TRUE) > 50 || min(fpkm_data, na.rm = TRUE) < 0) {
  cat("Applying log2(FPKM+1) transformation\n")
  expr_log <- log2(fpkm_data + 1)
} else {
  cat("Assuming data is already log-transformed\n")
  expr_log <- fpkm_data
}

# -------------------------- Create SingleCellExperiment object --------------------------
sce <- SingleCellExperiment(assays = list(logcounts = expr_log), colData = cell_meta)

# -------------------------- PCA (required for UMAP initialization) --------------------------
log_mat_t <- t(expr_log)
pca_res <- prcomp(log_mat_t, scale. = TRUE, center = TRUE)
reducedDim(sce, "PCA") <- pca_res$x

# -------------------------- UMAP --------------------------
n_neighbors <- min(20, n_samples - 1)
cat("UMAP using n_neighbors =", n_neighbors, "\n")
sce <- runUMAP(sce, dimred = "PCA", n_neighbors = n_neighbors)

# -------------------------- Extract UMAP coordinates --------------------------
umap_coords <- reducedDim(sce, "UMAP")
colnames(umap_coords) <- c("UMAP1", "UMAP2")
umap_df <- data.frame(sample = sample_names, group = cell_meta$group, umap_coords)

# -------------------------- Save UMAP coordinates --------------------------
write.xlsx(umap_df, file.path(output_dir, "umap_coords.xlsx"), rowNames = FALSE)

# -------------------------- Visualization --------------------------
unique_groups <- unique(cell_meta$group)
n_groups <- length(unique_groups)

safe_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", 
                 "#FF7F00", "#FFFF33", "#A65628", "#F781BF")
if (n_groups <= length(safe_colors)) {
  group_colors <- setNames(safe_colors[1:n_groups], unique_groups)
} else {
  group_colors <- setNames(rainbow(n_groups), unique_groups)
}

p_umap_label <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = group, label = sample)) +
  geom_point(size = 3) + scale_color_manual(values = group_colors) +
  geom_text_repel(size = 3, max.overlaps = Inf) +
  labs(title = "UMAP (FPKM)") + theme_minimal()
ggsave(file.path(output_dir, "UMAP_with_labels.pdf"), p_umap_label, width = 8, height = 6)

p_umap <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = group)) +
  geom_point(size = 3) + scale_color_manual(values = group_colors) +
  labs(title = "UMAP (FPKM)") + theme_minimal()
ggsave(file.path(output_dir, "UMAP.pdf"), p_umap, width = 7, height = 5)

cat("\nAnalysis complete! UMAP results saved in:", output_dir, "\n")
