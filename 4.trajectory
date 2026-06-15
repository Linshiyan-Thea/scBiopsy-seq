# ==============================================================================
# Pseudotime Trajectory Analysis via Principal Curve
# Description: Perform dimensionality reduction (PCA, UMAP, t-SNE) and fit a
#              principal curve to infer pseudotemporal ordering of single-cell
#              samples across two time points.
# Input:  - A gene expression matrix (genes x samples) in Excel format.
#         - A sample metadata file with sample names and group labels.
# Output: Trajectory plots (UMAP & t-SNE), pseudotime tables, and gene
#         filtering summaries saved to the specified output directory.
# ==============================================================================

library(openxlsx)
library(dplyr)
library(stringr)
library(SingleCellExperiment)
library(scater)
library(ggplot2)
library(ggrepel)
library(viridis)
library(princurve)

# ======================== User settings ======================================
# Replace with your own file paths
count_file <- "path/to/your/count_matrix.xlsx"
group_file <- "path/to/your/sample_groups.xlsx"
output_dir <- "path/to/your/output"

# Gene filtering: retain genes expressed in at least this percentage of samples
filter_percent <- 70

# Specify start and end groups to orient pseudotime direction
start_group <- "T0"
end_group   <- "T1"

# Whether to flip pseudotime so that start_group has the lowest values
adjust_direction <- TRUE

# ======================== Load expression data ================================
cat("Loading expression matrix...\n")
df_raw <- read.xlsx(count_file, rowNames = FALSE, check.names = FALSE)

gene_names <- as.character(df_raw[, 1])
dup_flag   <- duplicated(gene_names)
if (any(dup_flag)) {
  cat("Found", sum(dup_flag), "duplicate gene names; keeping first occurrence only.\n")
  df_raw     <- df_raw[!dup_flag, ]
  gene_names <- gene_names[!dup_flag]
}
rownames(df_raw)   <- gene_names
count_data_raw     <- as.matrix(df_raw[, -1, drop = FALSE])
count_data_raw[is.na(count_data_raw)] <- 0

# Sanitize sample names
sample_names <- trimws(colnames(count_data_raw))
sample_names <- gsub("-", ".", sample_names)
colnames(count_data_raw) <- sample_names

# ======================== Load sample metadata ================================
cat("Loading sample metadata...\n")
meta_raw  <- read.xlsx(group_file)
col_names <- colnames(meta_raw)

sample_col <- ifelse("sample"  %in% col_names, "sample",
              ifelse("Sample"  %in% col_names, "Sample", col_names[1]))
group_col  <- ifelse("group"   %in% col_names, "group",
              ifelse("Group"   %in% col_names, "Group",  col_names[2]))

group_df <- data.frame(
  sample = trimws(gsub("-", ".", as.character(meta_raw[[sample_col]]))),
  group  = as.character(meta_raw[[group_col]]),
  stringsAsFactors = FALSE
)

if (!all(sample_names %in% group_df$sample)) {
  missing <- setdiff(sample_names, group_df$sample)
  stop("The following samples are missing from the metadata file: ",
       paste(missing, collapse = ", "))
}
group_df  <- group_df[match(sample_names, group_df$sample), ]
cell_meta <- data.frame(sample = sample_names, group = group_df$group,
                        row.names = sample_names)
cat("Group summary:\n")
print(table(cell_meta$group))

if (!start_group %in% cell_meta$group) stop("Start group not found in data.")
if (!end_group   %in% cell_meta$group) stop("End group not found in data.")

n_total_samples <- ncol(count_data_raw)

# ======================== Gene filtering =====================================
cat("\n========== Filter threshold: ", filter_percent, "% ==========\n", sep = "")

percent_dir <- file.path(output_dir, paste0("percent_", filter_percent))
if (!dir.exists(percent_dir)) dir.create(percent_dir, recursive = TRUE)

min_samples <- ceiling(filter_percent / 100 * n_total_samples)
keep_genes  <- apply(count_data_raw, 1, function(x) sum(x > 0) >= min_samples)
count_data  <- count_data_raw[keep_genes, , drop = FALSE]

cat("Minimum samples required: ", min_samples, "\n", sep = "")
cat("Genes before filtering: ", nrow(count_data_raw), "\n")
cat("Genes after  filtering: ", sum(keep_genes), "\n")

if (sum(keep_genes) == 0) stop("No genes retained. Lower the filter threshold.")

gene_summary <- data.frame(percent = filter_percent, n_genes = sum(keep_genes))
write.csv(gene_summary, file.path(percent_dir, "gene_count_summary.csv"),
          row.names = FALSE)

# ======================== Build SCE & normalize ==============================
sce <- SingleCellExperiment(assays = list(counts = count_data),
                            colData = cell_meta)
sce <- logNormCounts(sce)

# ======================== Dimensionality reduction ============================
log_mat   <- logcounts(sce)
pca_res   <- prcomp(t(log_mat), scale. = TRUE, center = TRUE)
reducedDim(sce, "PCA") <- pca_res$x

n_samples   <- ncol(count_data)
n_neighbors <- min(8, max(2, n_samples - 1))
cat("UMAP n_neighbors =", n_neighbors, "\n")
sce <- runUMAP(sce, dimred = "PCA", n_neighbors = n_neighbors)

perplexity <- max(6, min(30, floor((n_samples - 1) / 3)))
cat("t-SNE perplexity =", perplexity, "\n")
sce <- runTSNE(sce, dimred = "PCA", perplexity = perplexity)

# ======================== Principal curve fitting =============================
umap_coords <- reducedDim(sce, "UMAP"); colnames(umap_coords) <- c("UMAP1", "UMAP2")
fit_umap    <- principal_curve(umap_coords, smoother = "smooth.spline", maxit = 50)
pseudotime_umap <- fit_umap$lambda

tsne_coords <- reducedDim(sce, "TSNE"); colnames(tsne_coords) <- c("tSNE1", "tSNE2")
fit_tsne    <- principal_curve(tsne_coords, smoother = "smooth.spline", maxit = 50)
pseudotime_tsne <- fit_tsne$lambda

# ======================== Orient pseudotime direction =========================
if (adjust_direction) {
  group_vec <- cell_meta$group

  # UMAP
  mean_start <- mean(pseudotime_umap[group_vec == start_group], na.rm = TRUE)
  mean_end   <- mean(pseudotime_umap[group_vec == end_group],   na.rm = TRUE)
  if (mean_start > mean_end) {
    pseudotime_umap <- max(pseudotime_umap) - pseudotime_umap
    cat("UMAP pseudotime reversed (start group mean > end group mean).\n")
  }

  # t-SNE
  mean_start <- mean(pseudotime_tsne[group_vec == start_group], na.rm = TRUE)
  mean_end   <- mean(pseudotime_tsne[group_vec == end_group],   na.rm = TRUE)
  if (mean_start > mean_end) {
    pseudotime_tsne <- max(pseudotime_tsne) - pseudotime_tsne
    cat("t-SNE pseudotime reversed.\n")
  }

  # Scale to [0, 1]
  pseudotime_umap <- (pseudotime_umap - min(pseudotime_umap)) /
                     (max(pseudotime_umap) - min(pseudotime_umap))
  pseudotime_tsne <- (pseudotime_tsne - min(pseudotime_tsne)) /
                     (max(pseudotime_tsne) - min(pseudotime_tsne))
}

# ======================== Assemble output data frames =========================
umap_df <- data.frame(
  sample     = sample_names,
  group      = cell_meta$group,
  pseudotime = pseudotime_umap,
  UMAP1      = umap_coords[, 1],
  UMAP2      = umap_coords[, 2]
)

tsne_df <- data.frame(
  sample     = sample_names,
  group      = cell_meta$group,
  pseudotime = pseudotime_tsne,
  tSNE1      = tsne_coords[, 1],
  tSNE2      = tsne_coords[, 2]
)

write.xlsx(umap_df, file.path(percent_dir, "trajectory_umap.xlsx"), rowNames = FALSE)
write.xlsx(tsne_df, file.path(percent_dir, "trajectory_tsne.xlsx"), rowNames = FALSE)

# ======================== Ordered curve points (for arrows) ===================
umap_curve_df <- data.frame(
  UMAP1 = fit_umap$s[order(pseudotime_umap), 1],
  UMAP2 = fit_umap$s[order(pseudotime_umap), 2]
)

tsne_curve_df <- data.frame(
  tSNE1 = fit_tsne$s[order(pseudotime_tsne), 1],
  tSNE2 = fit_tsne$s[order(pseudotime_tsne), 2]
)

# ======================== Visualization =======================================
arrow_style <- arrow(angle = 20, length = unit(0.3, "inches"), type = "closed")

# Define group colors — adjust names to match your group labels
group_colors <- c("T0" = "blue", "T1" = "red")

plot_trajectory <- function(df, x_var, y_var, curve_df, cx_var, cy_var,
                            color_var, palette, title, filename,
                            add_labels = FALSE, dashed = FALSE) {
  p <- ggplot(df, aes_string(x = x_var, y = y_var, color = color_var)) +
    geom_point(size = 3) +
    geom_path(data = curve_df, aes_string(x = cx_var, y = cy_var),
              color = "black", linewidth = 0.8,
              linetype = ifelse(dashed, "dashed", "solid"),
              arrow = arrow_style)

  if (color_var == "pseudotime") {
    p <- p + scale_color_viridis_c(option = "plasma")
  } else {
    p <- p + scale_color_manual(values = palette, name = "Group")
  }

  if (add_labels) {
    p <- p + geom_text_repel(aes(label = sample), size = 3, max.overlaps = Inf,
                             box.padding = 0.5, point.padding = 0.3)
  }

  p <- p + labs(title = title, color = ifelse(color_var == "pseudotime",
                                               "Pseudotime", "Group")) +
    theme_minimal()

  ggsave(file.path(percent_dir, filename), p, width = 8, height = 6)
}

base_title_umap <- paste0("Trajectory (UMAP) - ", filter_percent, "% expressed genes")
base_title_tsne <- paste0("Trajectory (t-SNE) - ", filter_percent, "% expressed genes")

# UMAP plots
plot_trajectory(umap_df, "UMAP1", "UMAP2", umap_curve_df, "UMAP1", "UMAP2",
                "pseudotime", NULL, paste0("Pseudotime ", base_title_umap),
                "trajectory_umap_pseudotime_with_labels.pdf", add_labels = TRUE)
plot_trajectory(umap_df, "UMAP1", "UMAP2", umap_curve_df, "UMAP1", "UMAP2",
                "pseudotime", NULL, paste0("Pseudotime ", base_title_umap),
                "trajectory_umap_pseudotime.pdf")
plot_trajectory(umap_df, "UMAP1", "UMAP2", umap_curve_df, "UMAP1", "UMAP2",
                "group", group_colors, base_title_umap,
                "trajectory_umap_group_with_labels.pdf", add_labels = TRUE)
plot_trajectory(umap_df, "UMAP1", "UMAP2", umap_curve_df, "UMAP1", "UMAP2",
                "group", group_colors, base_title_umap,
                "trajectory_umap_group.pdf")

# t-SNE plots
plot_trajectory(tsne_df, "tSNE1", "tSNE2", tsne_curve_df, "tSNE1", "tSNE2",
                "pseudotime", NULL, paste0("Pseudotime ", base_title_tsne),
                "trajectory_tsne_pseudotime_with_labels.pdf",
                add_labels = TRUE, dashed = TRUE)
plot_trajectory(tsne_df, "tSNE1", "tSNE2", tsne_curve_df, "tSNE1", "tSNE2",
                "pseudotime", NULL, paste0("Pseudotime ", base_title_tsne),
                "trajectory_tsne_pseudotime.pdf", dashed = TRUE)
plot_trajectory(tsne_df, "tSNE1", "tSNE2", tsne_curve_df, "tSNE1", "tSNE2",
                "group", group_colors, base_title_tsne,
                "trajectory_tsne_group_with_labels.pdf",
                add_labels = TRUE, dashed = TRUE)
plot_trajectory(tsne_df, "tSNE1", "tSNE2", tsne_curve_df, "tSNE1", "tSNE2",
                "group", group_colors, base_title_tsne,
                "trajectory_tsne_group.pdf", dashed = TRUE)

# ======================== Save pseudotime tables ==============================
write.csv(umap_df[, c("sample", "group", "pseudotime")],
          file.path(percent_dir, "pseudotime_umap.csv"), row.names = FALSE)
write.csv(tsne_df[, c("sample", "group", "pseudotime")],
          file.path(percent_dir, "pseudotime_tsne.csv"), row.names = FALSE)

cat("\nDone. Results saved to: ", percent_dir, "\n")
