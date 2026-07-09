# ==============================================================================
# Correlation Analysis Between Two Methods
# Description: Compute Pearson and Spearman correlations of mean gene expression
#              between two experimental conditions and visualise as scatter plots.
# Input:  - A gene expression matrix (genes x samples) in Excel format.
#         - A sample metadata file with columns: sample, group.
# Output: Scatter plots (PDF) and correlation coefficients (txt).
# ==============================================================================
# Environment: R 4.4.3
# Packages: ggplot2 4.0.2, readxl 1.4.5
# ==============================================================================

library(ggplot2)
library(readxl)

# ======================== User settings ======================================
# Replace with your own file paths
matrix_path <- "path/to/your/expression_matrix.xlsx"
class_path  <- "path/to/your/sample_metadata.xlsx"
out_dir     <- "path/to/your/output"

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ======================== Load data ==========================================
# Expression matrix: first column = Gene, remaining columns = samples
expr_df <- read_excel(matrix_path)
rownames(expr_df) <- expr_df[[1]]
expr_df <- expr_df[, -1]

# Sample metadata: columns assumed to be sample, group
class_info <- read_excel(class_path)
colnames(class_info) <- c("sample", "group")

# ======================== Match samples & split groups =======================
samples <- colnames(expr_df)

if (!all(samples %in% class_info$sample)) {
  missing <- samples[!samples %in% class_info$sample]
  warning("Samples not found in metadata: ", paste(missing, collapse = ", "))
}

group_vec <- class_info$group[match(samples, class_info$sample)]

# Auto-detect the two groups
group_levels <- unique(na.omit(group_vec))
if (length(group_levels) != 2) {
  stop("The 'group' column must contain exactly two levels. Found: ",
       paste(group_levels, collapse = ", "))
}
group_a <- group_levels[1]
group_b <- group_levels[2]
cat("Groups detected: ", group_a, " vs ", group_b, "\n")

idx_a <- which(group_vec == group_a)
idx_b <- which(group_vec == group_b)

if (length(idx_a) == 0) stop("No samples found for group: ", group_a)
if (length(idx_b) == 0) stop("No samples found for group: ", group_b)

mat_a <- expr_df[, idx_a, drop = FALSE]
mat_b <- expr_df[, idx_b, drop = FALSE]

# ======================== Compute mean & log-transform =======================
mean_a <- rowMeans(mat_a, na.rm = TRUE)
mean_b <- rowMeans(mat_b, na.rm = TRUE)

log2_a <- log2(mean_a + 0.25)
log2_b <- log2(mean_b + 0.25)

plot_data <- na.omit(data.frame(x = log2_a, y = log2_b))

# ======================== Correlation analysis ================================
methods <- c("pearson", "spearman")

for (method in methods) {
  cor_value <- cor(plot_data$x, plot_data$y, method = method)
  gene_num  <- nrow(plot_data)

  result_df <- data.frame(method = method, correlation = cor_value, n_genes = gene_num)
  write.table(result_df,
              file = file.path(out_dir, paste0("correlation_", method, ".txt")),
              col.names = TRUE, row.names = FALSE, sep = "\t")

  cat(toupper(method), " R = ", round(cor_value, 3), "\n")
  cat("Genes used = ", gene_num, "\n\n")

  # ======================== Scatter plot =======================================
  p <- ggplot(plot_data, aes(x = x, y = y)) +
    geom_point(alpha = 0.6, size = 1.5) +
    geom_smooth(method = "lm", se = FALSE, color = "red") +
    geom_text(
      label = paste0(ifelse(method == "pearson", "R = ", "rho = "),
                     round(cor_value, 3)),
      x = min(plot_data$x) + 0.05 * diff(range(plot_data$x)),
      y = min(plot_data$y) + 0.92 * diff(range(plot_data$y)),
      size = 5
    ) +
    labs(
      title = paste0(toupper(method), " correlation"),
      x = paste0("Log2(expression + 0.25) of ", group_a),
      y = paste0("Log2(expression + 0.25) of ", group_b)
    ) +
    theme_bw() +
    theme(
      panel.border     = element_blank(),
      axis.line        = element_line(color = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title       = element_text(hjust = 0.5)
    )

  pdf(file = file.path(out_dir, paste0("correlation_", method, ".pdf")),
      width = 4, height = 4)
  print(p)
  dev.off()

  print(p)
}

cat("\nDone. Results saved to: ", out_dir, "\n")
