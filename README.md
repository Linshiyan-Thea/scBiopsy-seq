# scBiopsy-seq

scBiopsy-seq: a temporal scRNA-seq assay combining electroosmosis-electrophoresis extraction & digital microfluidics. Detects 10–15K genes/cell with >90% success, enables sequential cytoplasmic extraction to link phenotype with transcriptional dynamics in development, viral infection & transcriptional suppression.

## Overview

This repository contains the computational pipeline for scBiopsy-seq data analysis. The workflow is organized into six sequential steps:

1. **Preprocessing** — quality control, adapter trimming, alignment, and read counting.
2. **Clustering** — UMAP clustering based on exon expression, within-cell expression changes, and cross-platform comparison.
3. **Score calculation** — geometric-mean-based scoring for gene sets of interest.
4. **Trajectory analysis** — pseudotime inference via Slingshot with gene–pseudotime correlation.
5. **Differential expression** — DESeq2-based DE analysis supporting paired and unpaired designs.
6. **Correlation analysis** — Pearson/Spearman correlation of mean gene expression between two conditions.


## Repository structure

```
scBiopsy-seq/
├── README.md
├── 1.preprocessing/
│   ├── S0.raw_fastqc.py              # FastQC quality control of raw reads
│   ├── S1.adapters_RNA.fasta         # Adapter sequences for trimming
│   ├── S1.trim_fastqc_ver0.02.py     # Adapter trimming and post-trim QC
│   ├── S2.star_map_ver0.02.py        # STAR alignment
│   ├── S3.inf_star_map_ver0.1.py     # Aggregate STAR mapping statistics
│   ├── S4.ExonIntron_RNA_v2.2.py     # Exon/intron mapping rate calculation
│   ├── S5.inf_ExonIntron_v2.py       # Aggregate exon/intron statistics
│   └── S6.DepthGeneReadCount_v5.2a.py # Read-depth subsampling and gene counting
├── 2.clustering/
│   ├── cluster based on exon.R                  # UMAP clustering on exon FPKM
│   ├── clustering based on expression change.R  # UMAP + k-means on log2 FC
│   └── clustering based on scBiopsy and scRNA-seq.R  # Cross-platform UMAP
├── 3.score/
│   └── score_calculation.R           # Gene set score (geometric mean)
├── 4.trajectory/
│   └── Trajectory.R                  # Slingshot pseudotime trajectory
├── 5.DE_analysis/
│   └── DE_analysis.R                 # DESeq2 differential expression
└── 6.correlation/
    └── correlation.R                 # Inter-condition correlation
```


## Environment

### Python (preprocessing)

```
Python 3.12
pandas 3.0.3
numpy 2.5.0
openpyxl 3.1.5
```

External command-line tools:

| Tool | Version |
|------|---------|
| FastQC | v0.12.1 |
| cutadapt | 5.0 |
| STAR | 2.7.11b |
| samtools | 1.19.2 |

### R (downstream analysis)

```
R 4.4.3
```

Key packages:

| Package | Version | Used in |
|---------|---------|---------|
| DESeq2 | 1.46.0 | 5.DE_analysis |
| SingleCellExperiment | 1.28.1 | 2.clustering, 4.trajectory |
| scater | 1.34.1 | 2.clustering, 4.trajectory |
| slingshot | 2.14.0 | 4.trajectory |
| umap | 0.2.10.0 | 2.clustering |
| ggplot2 | 4.0.2 | 2–6 |
| ggrepel | 0.9.6 | 2, 4, 5 |
| openxlsx | 4.2.8 | 2, 3, 4, 5 |
| readxl | 1.4.5 | 3, 6 |
| cowplot | 1.2.0 | 2.clustering |
| tidyverse | 2.0.0 | 3.score |
| viridis | 0.6.5 | 4.trajectory |
| patchwork | 1.3.2 | 4.trajectory |
| dplyr | 1.1.4 | 2, 4, 5 |
| stringr | 1.5.2 | 4.trajectory |


## Usage

### 1. Preprocessing

Run scripts sequentially in the `1.preprocessing/` directory. Each script accepts command-line arguments; use `--help` for details.

```bash
# Quality control
python S0.raw_fastqc.py

# Trimming
python S1.trim_fastqc_ver0.02.py

# STAR alignment
python S2.star_map_ver0.02.py --star_index path/to/STAR/genome/index --gtf path/to/annotation.gtf

# Mapping statistics
python S3.inf_star_map_ver0.1.py -i ./log_final_out -o output_prefix

# Exon/intron rate
python S4.ExonIntron_RNA_v2.2.py -e path/to/exon.bed -g path/to/gene.bed

# Aggregate exon/intron info
python S5.inf_ExonIntron_v2.py -o output_prefix

# Subsampling and gene counting
python S6.DepthGeneReadCount_v5.2a.py -s list.sam -g path/to/annotation.gtf -mt path/to/mt_gene.txt -rp path/to/rp_gene_id.txt -d 1M -o output_prefix
```

Reference genome: **GRCh38.105** (human) or **GRCm38.102** (mouse).

### 2. Clustering

Edit the file paths at the top of each R script before running:

```r
# UMAP clustering based on exon expression
source("2.clustering/cluster based on exon.R")

# Clustering based on within-cell expression changes (log2 FC)
source("2.clustering/clustering based on expression change.R")

# Cross-platform comparison (scBiopsy-seq vs scRNA-seq)
source("2.clustering/clustering based on scBiopsy and scRNA-seq.R")
```

### 3. Gene set score

Provide a gene expression matrix containing only the genes of interest:

```r
source("3.score/score_calculation.R")
```

Set `gene_set_name` in the script to label the output file.

### 4. Trajectory analysis

```r
source("4.trajectory/Trajectory.R")
```

Set `start_group` to define the trajectory origin. Outputs include pseudotime tables, gene–pseudotime correlations, and top-10 gene UMAP plots.

### 5. Differential expression

```r
source("5.DE_analysis/DE_analysis.R")
```

Set `paired = TRUE/FALSE` depending on experimental design. For paired analysis, ensure the metadata file contains the pairing column (default: `pair`).

### 6. Correlation

```r
source("6.correlation/correlation.R")
```

Computes Pearson and Spearman correlations of mean gene expression between two groups.


## License

MIT License
# scBiopsy-seq
scBiopsy-seq: a temporal scRNA-seq assay combining electroosmosis-electrophoresis extraction &amp; digital microfluidics. Detects 10-15K genes/cell with >90% success, enables sequential cytoplasmic extraction to link phenotype with transcriptional dynamics in development, viral infection &amp; transcriptional suppression.
