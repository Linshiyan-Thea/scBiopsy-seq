# ==============================================================================
# IFN Response Score Calculation — Geometric Mean Method
# Description: Compute IFN response scores using the geometric mean of all
#              gene expression values across samples.
# Input:  An Excel file containing a gene expression matrix (genes x samples).
# Output: A CSV file with IFN scores for each sample.
# ==============================================================================

# 1. Load required packages ---------------------------------------------------
library(readxl)
library(tidyverse)

# 2. Read data ----------------------------------------------------------------
# Replace with your own file path
file_path <- "path/to/your/gene_expression.xlsx"
data_raw  <- read_excel(file_path, sheet = 1)

gene_names <- data_raw[[1]]
expr_mat   <- as.matrix(data_raw[, -1])
rownames(expr_mat) <- gene_names
colnames(expr_mat) <- colnames(data_raw)[-1]

# 3. Define groups ------------------------------------------------------------
# Adjust the number of samples per group according to your experimental design
n_control   <- 13
n_treatment <- 13
group <- c(rep("Control", n_control), rep("Treatment", n_treatment))

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
  Sample          = colnames(expr_mat),
  Group           = group,
  IFN_Score_geomean = geomean_scores
) %>% arrange(desc(IFN_Score_geomean))

cat("\n=== Geometric Mean Score ===\n")
print(result_geo)

# ==============================================================================
# Save results
# ==============================================================================
output_dir <- "path/to/your/output"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

write.csv(result_geo,
          file.path(output_dir, "IFN_scores_geomean.csv"),
          row.names = FALSE)

cat("\nResults saved to: ", output_dir, "\n")
