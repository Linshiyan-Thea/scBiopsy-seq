# ==============================================================================
# Gene Set Score Calculation — Geometric Mean Method
# Description: Compute a score for a given gene set using the geometric mean
#              of expression values across samples.
# Input:  - gene_expression.xlsx: gene expression matrix (genes x samples,
#           first column = gene names).
#         - sample_groups.xlsx: sample grouping file (columns: sample, group).
# Output: A CSV file with scores for each sample.
# ==============================================================================
# Environment: R 4.4.3
# Packages: readxl 1.4.5, tidyverse 2.0.0
# ==============================================================================

# 1. Load required packages ---------------------------------------------------
library(readxl)
library(tidyverse)

# -------------------------- User settings ------------------------------------
gene_set_name <- "my_gene_set"  # label for the gene set (used in output file names)

file_path  <- "path/to/your/gene_expression.xlsx"
group_file <- "path/to/your/sample_groups.xlsx"
output_dir <- "path/to/your/output"
# -----------------------------------------------------------------------------

# 2. Read data ----------------------------------------------------------------
data_raw  <- read_excel(file_path, sheet = 1)

gene_names <- data_raw[[1]]
expr_mat   <- as.matrix(data_raw[, -1])
rownames(expr_mat) <- gene_names
colnames(expr_mat) <- colnames(data_raw)[-1]

# 3. Read grouping information ------------------------------------------------
group_df <- read_excel(group_file)
colnames(group_df) <- c("sample", "group")
group_df$sample <- trimws(group_df$sample)

# Match sample order between expression matrix and grouping file
sample_names <- colnames(expr_mat)
if (!all(sample_names %in% group_df$sample)) {
  missing <- sample_names[!sample_names %in% group_df$sample]
  stop("Samples missing in grouping file: ", paste(missing, collapse = ", "))
}
group_df <- group_df[match(sample_names, group_df$sample), ]
group <- group_df$group

# ==============================================================================
# Geometric mean score
# Formula: exp(mean(log(x + pseudo))) / 10
# ==============================================================================
pseudo <- 0.01
expr_pseudo <- expr_mat + pseudo

geomean_scores <- apply(expr_pseudo, 2, function(x) {
  exp(mean(log(x))) / 10
})

result_geo <- data.frame(
  Sample = colnames(expr_mat),
  Group  = group,
  Score  = geomean_scores
) %>% arrange(desc(Score))

cat("\n=== Geometric Mean Score ===\n")
print(result_geo)

# ==============================================================================
# Save results
# ==============================================================================
output_dir <- "path/to/your/output"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

write.csv(result_geo,
          file.path(output_dir, paste0(gene_set_name, "_scores_geomean.csv")),
          row.names = FALSE)

cat("\nResults saved to: ", output_dir, "\n")
