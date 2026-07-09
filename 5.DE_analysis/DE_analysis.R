# ==============================================================================
# Differential Expression Analysis (DESeq2)
# Description: Perform DE analysis between two groups using DESeq2.
#              Supports both paired and unpaired designs.
# Input:  - A gene expression matrix (genes x samples) in Excel format.
#         - A sample metadata file with columns: sample, group (and optionally
#           a pairing column for paired designs).
# Output: DE results table and volcano plot.
# ==============================================================================
# Environment: R 4.4.3
# Packages: DESeq2 1.46.0, ggplot2 4.0.2, ggrepel 0.9.6,
#           openxlsx 4.2.8, dplyr 1.1.4
# ==============================================================================

# ======================== Parameters =========================================
padj_cutoff   <- 0.05   # Adjusted p-value threshold
log2FC_cutoff <- 1.5    # |log2 Fold Change| threshold
paired        <- TRUE   # Set to TRUE for paired analysis, FALSE for unpaired
pair_id_col   <- "pair" # Column name in colData for pairing (used if paired = TRUE)

# Replace with your own file paths
expr_file  <- "path/to/your/expression_matrix.xlsx"
meta_file  <- "path/to/your/sample_metadata.xlsx"
output_dir <- "path/to/your/output"

# ======================== Load libraries =====================================
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(openxlsx)
library(dplyr)


# ======================== Load expression matrix ==============================
cat("Loading expression matrix...\n")
count_data <- read.xlsx(expr_file, sheet = 1, rowNames = FALSE)
colnames(count_data)[1] <- "Gene"

# Remove empty gene names
count_data <- count_data %>% filter(!is.na(Gene) & Gene != "")

# Merge duplicate genes by summing
count_matrix <- count_data %>%
  group_by(Gene) %>%
  summarise(across(everything(), ~ sum(.x, na.rm = TRUE))) %>%
  ungroup() %>%
  as.data.frame()

rownames(count_matrix) <- count_matrix$Gene
count_matrix <- count_matrix[, -1]

# ======================== Load sample metadata ================================
cat("Loading sample metadata...\n")
colData <- read.xlsx(meta_file, sheet = 1)

# Sanitize sample names
colnames(count_matrix) <- tolower(make.names(colnames(count_matrix)))
colData$sample <- tolower(make.names(colData$sample))

# Align sample order
count_matrix <- count_matrix[, colData$sample]

# Auto-detect the two groups from colData$group
group_levels <- unique(colData$group)
if (length(group_levels) != 2) {
  stop("The 'group' column must contain exactly two levels. Found: ",
       paste(group_levels, collapse = ", "))
}
group_ref <- group_levels[1]   # Reference group (denominator)
group_alt <- group_levels[2]   # Comparison group (numerator)
cat("Groups detected: ", group_ref, " (reference) vs ", group_alt, "\n")

# ======================== Validate paired design =============================
if (paired) {
  if (!pair_id_col %in% colnames(colData)) {
    stop("Paired analysis requested (paired = TRUE), but column '", pair_id_col,
         "' not found in metadata. Available columns: ",
         paste(colnames(colData), collapse = ", "),
         "\nSet paired = FALSE or add the pairing column to your metadata file.")
  }
  cat("Paired design: using '", pair_id_col, "' as pairing variable.\n")
}

# ======================== Build DESeq2 object =================================
design_formula <- if (paired) {
  as.formula(paste0("~ ", pair_id_col, " + group"))
} else {
  ~ group
}

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData   = colData,
  design    = design_formula
)

cat("Design formula: ", deparse(design_formula), "\n")

# ======================== DE analysis =========================================
dds <- DESeq(dds)

res <- results(dds,
               contrast = c("group", group_alt, group_ref),
               alpha    = padj_cutoff)

DEG <- as.data.frame(res)
DEG$gene <- rownames(DEG)

# Label genes as Up / Down / Not significant
DEG$change <- as.factor(
  ifelse(DEG$padj < padj_cutoff & abs(DEG$log2FoldChange) >= log2FC_cutoff,
         ifelse(DEG$log2FoldChange > log2FC_cutoff, "Up", "Down"),
         "Not")
)

cat("\nDEG summary:\n")
print(table(DEG$change))

# ======================== Volcano plot ========================================
top_genes <- DEG[order(DEG$pvalue, -abs(DEG$log2FoldChange)), ][1:30, ]

x_abs_max <- max(abs(na.omit(DEG$log2FoldChange))) + 0.5
y_max     <- max(-log10(na.omit(DEG$padj)), na.rm = TRUE) + 0.5

volcano_plot <- ggplot(DEG, aes(x = log2FoldChange, y = -log10(padj),
                                color = change)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("Up"   = "#E95D49",
                                "Down" = "#5C63A4",
                                "Not"  = "grey")) +
  geom_text_repel(data = top_genes, aes(label = gene),
                  size = 3, max.overlaps = 20, box.padding = 0.5) +
  geom_hline(yintercept = -log10(padj_cutoff),
             linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-log2FC_cutoff, log2FC_cutoff),
             linetype = "dashed", color = "black") +
  labs(title = paste0("DEG: ", group_alt, " vs ", group_ref,
                      " (padj<", padj_cutoff,
                      ", |log2FC|>=", log2FC_cutoff, ")"),
       x = "log2 Fold Change",
       y = "-log10(adjusted p-value)") +
  scale_x_continuous(limits = c(-x_abs_max, x_abs_max), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, y_max), expand = c(0, 0)) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

print(volcano_plot)

# ======================== Save results ========================================
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

DEG_output <- cbind(gene = rownames(DEG), DEG)
write.csv(DEG_output,
          file = file.path(output_dir, "DESeq2_results.csv"),
          row.names = FALSE)

ggsave(filename = file.path(output_dir, "Volcano_plot.pdf"),
       plot = volcano_plot, width = 8, height = 6, dpi = 300)

cat("\nDone. Results saved to: ", output_dir, "\n")
