library(openxlsx)
library(ggplot2)
library(ggrepel)
library(umap)

count_file <- "Matrix.xlsx"
group_file <- "class_info_DESeq2.xlsx"
output_dir <- "UMAP_FPKM_results"

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Read expression matrix
df_raw <- read.xlsx(count_file, rowNames = FALSE, check.names = FALSE)
gene_names <- as.character(df_raw[, 1])

keep <- !is.na(gene_names)
df_raw <- df_raw[keep, ]
gene_names <- gene_names[keep]

dup_flag <- duplicated(gene_names)
if (any(dup_flag)) {
  df_raw <- df_raw[!dup_flag, ]
  gene_names <- gene_names[!dup_flag]
}

rownames(df_raw) <- gene_names
fpkm_data <- as.matrix(df_raw[, -1, drop = FALSE])
fpkm_data[is.na(fpkm_data)] <- 0

# Filter lowly expressed genes
keep_genes <- rowSums(fpkm_data > 0) >= min_samples
fpkm_data <- fpkm_data[keep_genes, ]

# Clean sample names
sample_names <- colnames(fpkm_data)
sample_names <- trimws(sample_names)
sample_names <- gsub("-", ".", sample_names)
colnames(fpkm_data) <- sample_names

# Read grouping
group_df <- read.xlsx(group_file)
colnames(group_df) <- c("sample", "group")
group_df$sample <- trimws(group_df$sample)
group_df$sample <- gsub("-", ".", group_df$sample)
group_df <- group_df[match(sample_names, group_df$sample), ]
cell_meta <- data.frame(sample = sample_names, group = group_df$group, stringsAsFactors = FALSE)

# Log transform if needed
if (max(fpkm_data, na.rm = TRUE) > 50) {
  expr_log <- log2(fpkm_data + 1)
} else {
  expr_log <- fpkm_data
}

# PCA then UMAP
pca_res <- prcomp(t(expr_log), scale. = TRUE, center = TRUE)
n_neighbors <- min(umap_neighbors, ncol(fpkm_data) - 1)
umap_res <- umap(pca_res$x[, 1:min(50, ncol(pca_res$x))],
                 n_neighbors = n_neighbors,
                 n_components = 2)
umap_df <- data.frame(sample = sample_names,
                      group = cell_meta$group,
                      UMAP1 = umap_res$layout[, 1],
                      UMAP2 = umap_res$layout[, 2])

# Save coordinates
write.xlsx(umap_df, file.path(output_dir, "umap_coords.xlsx"), rowNames = FALSE)

# Plot
group_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00")[1:length(unique(cell_meta$group))]
names(group_colors) <- unique(cell_meta$group)

p <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = group, label = sample)) +
  geom_point(size = 3) +
  scale_color_manual(values = group_colors) +
  geom_text_repel(size = 3, max.overlaps = Inf) +
  labs(title = "UMAP (FPKM)") +
  theme_minimal()
ggsave(file.path(output_dir, "UMAP_with_labels.pdf"), p, width = 8, height = 6)

p2 <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = group)) +
  geom_point(size = 3) +
  scale_color_manual(values = group_colors) +
  labs(title = "UMAP (FPKM)") +
  theme_minimal()
ggsave(file.path(output_dir, "UMAP.pdf"), p2, width = 7, height = 5)

cat("Done. Results in", output_dir, "\n")
