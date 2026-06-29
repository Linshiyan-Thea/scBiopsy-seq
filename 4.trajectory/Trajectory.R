# ==============================================================================
# Pseudotime Trajectory Analysis via Slingshot
# Description: Perform dimensionality reduction (PCA, UMAP) and infer
#              pseudotemporal ordering of samples using Slingshot.
#              Includes gene expression–pseudotime correlation analysis.
# Input:  - A gene expression matrix (genes x samples) in Excel format.
#         - A sample metadata file with sample names and group labels.
# Output: Trajectory plots (UMAP), pseudotime tables, gene–pseudotime
#         correlation results, and top-10 gene mapping plots.
# ==============================================================================

library(openxlsx)
library(dplyr)
library(stringr)
library(SingleCellExperiment)
library(scater)
library(slingshot)
library(ggplot2)
library(ggrepel)
library(viridis)
library(patchwork)

# ======================== User settings ======================================
# Replace with your own file paths
count_file <- "path/to/your/count_matrix.xlsx"
group_file <- "path/to/your/sample_groups.xlsx"
output_dir <- "path/to/your/output"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Specify the start group for trajectory inference
start_group <- "Initial"

# ======================== Load expression data ================================
cat("Loading expression matrix...\n")
df_raw <- read.xlsx(count_file, rowNames = TRUE, check.names = FALSE)
count_data <- as.matrix(df_raw)
count_data[is.na(count_data)] <- 0

# Remove genes with zero expression across all samples
keep_genes <- rowSums(count_data) > 0
count_data <- count_data[keep_genes, ]
cat("Remaining genes:", nrow(count_data), "\n")

# ======================== Load sample metadata ================================
cat("Loading sample metadata...\n")
group_info <- read.xlsx(group_file)
sample_names <- colnames(count_data)

if (!all(sample_names %in% group_info$sample)) {
  stop("Sample names in the count matrix do not fully match the metadata file.")
}
group_info <- group_info[match(sample_names, group_info$sample), , drop = FALSE]

cell_meta <- data.frame(
  sample = sample_names,
  group  = group_info$group,
  row.names = sample_names
)
cat("Group summary:\n")
print(table(cell_meta$group))

if (!start_group %in% cell_meta$group) stop("Start group not found in data.")

# ======================== Build SCE & normalize ===============================
sce <- SingleCellExperiment(assays = list(counts = count_data),
                            colData = cell_meta)
sce <- logNormCounts(sce)  # log2(CPM + 1)

# ======================== Dimensionality reduction ============================
log_mat   <- logcounts(sce)
pca_res   <- prcomp(t(log_mat), scale. = TRUE, center = TRUE)
reducedDim(sce, "PCA") <- pca_res$x

n_samples   <- ncol(count_data)
n_neighbors <- min(15, n_samples - 1)
cat("UMAP n_neighbors =", n_neighbors, "\n")
sce <- runUMAP(sce, dimred = "PCA", n_neighbors = n_neighbors)

# ======================== Slingshot trajectory inference ======================
colData(sce)$cluster <- as.factor(cell_meta$group)
sce <- slingshot(sce, clusterLabels = 'cluster', reducedDim = 'UMAP',
                 start.clus = start_group)

# Extract pseudotime (first lineage)
pseudotime_vals <- slingPseudotime(sce)[, 1]

umap_coords <- reducedDim(sce, "UMAP")
colnames(umap_coords) <- c("UMAP1", "UMAP2")

# Assemble output data frame
plot_data <- data.frame(
  sample     = sample_names,
  pseudotime = pseudotime_vals,
  group      = cell_meta$group,
  UMAP1      = umap_coords[, 1],
  UMAP2      = umap_coords[, 2]
)

# Extract trajectory curves (for drawing arrows)
curves <- slingCurves(sce)
if (length(curves) > 0) {
  curves_list <- lapply(seq_along(curves), function(i) {
    curve <- curves[[i]]
    df <- as.data.frame(curve$s[curve$ord, ])
    colnames(df) <- c("UMAP1", "UMAP2")
    df$curve_id <- paste0("Lineage", i)
    df
  })
  curves_df <- do.call(rbind, curves_list)
} else {
  curves_df <- NULL
}

# Save trajectory data
write.csv(plot_data, file.path(output_dir, "pseudotime_results.csv"),
          row.names = FALSE)
write.xlsx(plot_data, file.path(output_dir, "trajectory_data_with_umap.xlsx"),
           rowNames = FALSE)

# ======================== Gene–pseudotime correlation =========================
cat("Computing gene–pseudotime correlations...\n")
log_counts_mat <- as.matrix(logcounts(sce))

# Remove samples with NA pseudotime (if any)
valid_idx <- which(!is.na(pseudotime_vals))
pseudotime_valid <- pseudotime_vals[valid_idx]
log_counts_valid <- log_counts_mat[, valid_idx, drop = FALSE]

# Spearman correlation (vectorized)
rho   <- apply(log_counts_valid, 1, function(x) cor(x, pseudotime_valid, method = "spearman", use = "complete.obs"))
p_val <- apply(log_counts_valid, 1, function(x) cor.test(x, pseudotime_valid, method = "spearman", exact = FALSE)$p.value)
p_adj <- p.adjust(p_val, method = "fdr")

cor_results <- data.frame(
  Gene         = rownames(log_counts_valid),
  Spearman_rho = rho,
  P_value      = p_val,
  FDR          = p_adj,
  stringsAsFactors = FALSE
)
cor_results <- cor_results[order(-abs(cor_results$Spearman_rho)), ]

# Save all genes
write.xlsx(cor_results, file.path(output_dir, "gene_pseudotime_correlation_all.xlsx"),
           rowNames = FALSE)

# Significant genes (|rho| > 0.5 & FDR < 0.05)
sig_genes <- cor_results[abs(cor_results$Spearman_rho) > 0.5 & cor_results$FDR < 0.05, ]
write.xlsx(sig_genes, file.path(output_dir, "gene_pseudotime_correlation_significant.xlsx"),
           rowNames = FALSE)

# Top 10 genes
top10_genes <- cor_results[1:min(10, nrow(cor_results)), "Gene"]
cat("Top 10 genes by |rho|:\n")
print(top10_genes)

# ======================== Top-10 gene mapping plots ===========================
umap_plot_df <- data.frame(umap_coords, pseudotime = pseudotime_vals,
                           group = cell_meta$group, sample = sample_names)

expr_top10   <- as.data.frame(t(log_counts_mat[top10_genes, , drop = FALSE]))
expr_top10$sample <- rownames(expr_top10)
plot_df <- merge(umap_plot_df, expr_top10, by = "sample")

# Define group colors — adjust names to match your group labels
group_colors <- c("Initial" = "blue", "Infected" = "red")

# Multi-page PDF: each page shows one gene (UMAP expression + trend)
pdf(file.path(output_dir, "top10_genes_mapping_plots.pdf"), width = 12, height = 6)
for (gene in top10_genes) {
  rho_val <- cor_results$Spearman_rho[cor_results$Gene == gene]
  p_fmt   <- format(cor_results$P_value[cor_results$Gene == gene], scientific = TRUE, digits = 3)
  fdr_fmt <- format(cor_results$FDR[cor_results$Gene == gene], scientific = TRUE, digits = 3)

  # Sub-plot 1: UMAP colored by gene expression
  p1 <- ggplot(plot_df, aes(x = UMAP1, y = UMAP2, color = .data[[gene]])) +
    geom_point(size = 2.5) +
    scale_color_viridis_c(option = "magma", name = "log2(CPM+1)") +
    labs(title = paste0(gene, " (rho = ", round(rho_val, 3), ", p = ", p_fmt, ")")) +
    theme_minimal() + theme(legend.position = "right")

  # Add trajectory curve arrow if available
  if (!is.null(curves_df)) {
    p1 <- p1 + geom_path(data = curves_df, aes(x = UMAP1, y = UMAP2, group = curve_id),
                         color = "black", linewidth = 0.6,
                         arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
                         inherit.aes = FALSE)
  }

  # Sub-plot 2: Expression vs pseudotime with loess smooth
  p2 <- ggplot(plot_df, aes(x = pseudotime, y = .data[[gene]])) +
    geom_point(aes(color = group), size = 2, alpha = 0.7) +
    scale_color_manual(values = group_colors, name = "Group") +
    geom_smooth(method = "loess", se = TRUE, color = "black", fill = "grey80") +
    labs(x = "Pseudotime", y = "log2(CPM+1)",
         title = paste0("Expression trend along pseudotime\nrho = ", round(rho_val, 3), ", FDR = ", fdr_fmt)) +
    theme_minimal() + theme(legend.position = "bottom")

  # Combine side by side
  combined <- p1 + p2 + plot_annotation(title = gene, theme = theme(plot.title = element_text(hjust = 0.5)))
  print(combined)
}
dev.off()
cat("Top-10 gene mapping plots saved.\n")

# ======================== Trajectory visualization ============================
arrow_style <- arrow(angle = 20, length = unit(0.3, "inches"), type = "closed")

# Plot 1: Pseudotime-colored + sample labels + trajectory arrow
p_traj <- ggplot(plot_data, aes(x = UMAP1, y = UMAP2, color = pseudotime, label = sample)) +
  geom_point(size = 3) +
  scale_color_viridis_c(option = "plasma") +
  geom_text_repel(size = 3, max.overlaps = Inf, box.padding = 0.5, point.padding = 0.3) +
  labs(title = "Pseudotime trajectory (UMAP)", color = "Pseudotime") +
  theme_minimal()
if (!is.null(curves_df)) {
  p_traj <- p_traj + geom_path(
    data = curves_df, aes(x = UMAP1, y = UMAP2, group = curve_id),
    color = "black", linewidth = 0.8,
    arrow = arrow_style, inherit.aes = FALSE
  )
}
ggsave(file.path(output_dir, "trajectory_umap_pseudotime_with_labels.pdf"),
       p_traj, width = 8, height = 6)

# Plot 2: Group-colored + sample labels + trajectory arrow
p_group <- ggplot(plot_data, aes(x = UMAP1, y = UMAP2, color = group, label = sample)) +
  geom_point(size = 3) +
  scale_color_manual(values = group_colors, name = "Group") +
  geom_text_repel(size = 3, max.overlaps = Inf, box.padding = 0.5, point.padding = 0.3) +
  labs(title = "Trajectory colored by group (UMAP)") +
  theme_minimal()
if (!is.null(curves_df)) {
  p_group <- p_group + geom_path(
    data = curves_df, aes(x = UMAP1, y = UMAP2, group = curve_id),
    color = "black", linewidth = 0.8,
    arrow = arrow_style, inherit.aes = FALSE
  )
}
ggsave(file.path(output_dir, "trajectory_umap_by_group.pdf"),
       p_group, width = 8, height = 6)
