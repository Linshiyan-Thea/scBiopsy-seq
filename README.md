# scBiopsy-seq

scBiopsy-seq: a temporal scRNA-seq assay combining electroosmosis-electrophoresis extraction & digital microfluidics. Detects 10–15K genes/cell with ~96% success, enables sequential cytoplasmic extraction to link phenotype with transcriptional dynamics in development, viral infection & transcriptional suppression.

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
│   ├── S1.trim_fastqc.py             # Adapter trimming and post-trim QC
│   ├── S2.star_map.py                # STAR alignment
│   ├── S3.inf_star_map.py            # Aggregate STAR mapping statistics
│   ├── S4.ExonIntron_RNA.py          # Exon/intron mapping rate calculation
│   ├── S5.inf_ExonIntron.py          # Aggregate exon/intron statistics
│   └── S6.DepthGeneReadCount.py      # Read-depth subsampling and gene counting
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


## Installation

Clone the repository and enter the project directory:

```bash
git clone https://github.com/Linshiyan-Thea/scBiopsy-seq.git
cd scBiopsy-seq
```

### Python environment (preprocessing)

Create a conda environment for the preprocessing scripts:

```bash
conda create -n scbiopsy python=3.12 -y
conda activate scbiopsy
pip install pandas==3.0.3 numpy==2.5.0 openpyxl==3.1.5
```

The preprocessing pipeline also requires the following external command-line tools:

```bash
conda install -c bioconda fastqc=0.11.7 cutadapt=3.4 bbmap=38.90 star=2.7.3a samtools=1.3.1 htseq=0.12.4
```

> **Note:** The preprocessing scripts (S0–S6) invoke external tools via `os.system` and shell commands. These scripts must be run on a **Linux/Unix** system with the above tools available in `$PATH`.

### R environment (downstream analysis)

Install the R packages used for clustering, scoring, trajectory, DE, and correlation analysis:

```r
# CRAN packages
install.packages(c(
  "ggplot2", "ggrepel", "cowplot", "viridis", "patchwork",
  "dplyr", "stringr", "tidyr", "readxl", "openxlsx",
  "tidyverse", "umap"
))

# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c(
  "DESeq2", "SingleCellExperiment", "scater", "slingshot",
  "SC3", "clusterProfiler"
))
```

Package versions used in the published study are listed in each script's header comment.


## Environment

The versions listed below match those used in the published study.

### Python (preprocessing)

```
Python 3.12
pandas 3.0.3
numpy 2.5.0
openpyxl 3.1.5
```

### R (downstream analysis)

```
R 4.4.3
```

R package versions used in the published study are annotated in the header of each script.


## Usage

### 1. Preprocessing

> **Note:** All preprocessing scripts invoke external tools (FastQC, cutadapt, BBMap `repair.sh`, STAR, samtools, htseq-count) via `os.system` shell commands. They must be run on a **Linux/Unix** system with these tools available in `$PATH`.
>
> Scripts S0–S1 do not accept command-line arguments; edit the paths inside the script or place input files in the working directory as described below.
> Scripts S2–S6 accept command-line arguments; use `--help` for full parameter details.

#### S0 — Raw read quality control

Place all raw `*.fq` files in the working directory. The script runs FastQC on every `.fq` file and writes results to `fastqc_result/`.

```bash
python S0.raw_fastqc.py
```

#### S1 — Adapter trimming and post-trim QC

Place paired-end FASTQ files (`*1.fq`, `*2.fq`) and the adapter file `S1.adapters_RNA.fasta` in the working directory. The script runs cutadapt (adapter removal + poly-tail trimming), BBMap `repair.sh` (paired-end repair), and FastQC on trimmed reads.

```bash
python S1.trim_fastqc.py
```

Output: trimmed `*_repair_1.fq` / `*_repair_2.fq` files (used by S2).

#### S2 — STAR alignment

Align trimmed paired-end reads with STAR. Requires the STAR genome index and GTF annotation.

```bash
python S2.star_map.py --star_index path/to/STAR/genome/index --gtf path/to/annotation.gtf
```

Optional: `-t/--threads` (default: 10). Input: `*_1_repair_1.fq` / `*_2_repair_2.fq` from S1. Output: `Aligned.out.sam` and `Log.final.out` in `star_out_*/`.

#### S3 — Aggregate STAR mapping statistics

Collect all `*Log.final.out` files from the directory specified by `-i` and summarize mapping rates.

```bash
python S3.inf_star_map.py -i ./log_final_out -o output_prefix
```

`-i` defaults to `./log_final_out` if omitted.

#### S4 — Exon/intron mapping rate

Calculate exon and intron mapping rates for all `Aligned.out.sam` files in the working directory. Requires exon and gene BED files derived from the reference annotation.

```bash
python S4.ExonIntron_RNA.py -e path/to/exon.bed -g path/to/gene.bed
```

Output: `*_ExonIntron.txt` per sample (used by S5).

#### S5 — Aggregate exon/intron statistics

Collect all `*_ExonIntron.txt` files in the working directory and write a summary Excel table.

```bash
python S5.inf_ExonIntron.py -o output_prefix
```

Output: `output_prefix_ExonIntron_MAP.xlsx`.

#### S6 — Read-depth subsampling and gene counting

Subsample all samples to the same sequencing depth and count genes with htseq-count.

```bash
python S6.DepthGeneReadCount.py \
    -s sample1_Aligned.out.sam sample2_Aligned.out.sam \
    -g path/to/annotation.gtf \
    -mt path/to/mt_gene.txt \
    -rp path/to/rp_gene_id.txt \
    -d 1M \
    -o output_prefix
```

| Argument | Required | Description |
|----------|----------|-------------|
| `-s` | yes | One or more `.sam` files, or a text file listing sam paths |
| `-g` | yes | GTF annotation file |
| `-mt` | yes | Text file of mitochondrial gene IDs |
| `-rp` | yes | Text file of ribosomal protein gene IDs |
| `-d` | yes | Subsampling depth (e.g. `1M`) or `raw` to skip |
| `-o` | yes | Output prefix for `DepthGeneReadCount.xlsx` |
| `--seed` | no | Random seed for downsampling (default: 42) |

Reference genome: **GRCh38.105** (human) or **GRCm38.102** (mouse).

### 2. Clustering

Edit the file paths at the top of each R script before running. All scripts require the user to set `input_dir`, `output_dir`, and other parameters in the "User settings" block at the top.

```r
# UMAP clustering based on exon FPKM
source("2.clustering/cluster based on exon.R")

# UMAP + k-means clustering based on within-cell expression changes (log2 FC)
source("2.clustering/clustering based on expression change.R")

# Cross-platform comparison (scBiopsy-seq vs scRNA-seq)
source("2.clustering/clustering based on scBiopsy and scRNA-seq.R")
```

### 3. Gene set score

Compute a geometric-mean-based score for a custom gene set. Provide a gene expression matrix containing only the genes of interest and a sample group file.

```r
source("3.score/score_calculation.R")
```

Set `gene_set_name` in the script to label the output file. Set `input_file`, `group_file`, and `output_dir` in the "User settings" block.

### 4. Trajectory analysis

Perform Slingshot pseudotime trajectory inference with gene–pseudotime correlation analysis.

```r
source("4.trajectory/Trajectory.R")
```

Set `start_group` to define the trajectory origin. Edit `input_file`, `output_dir`, and `n_neighbors` in the "User settings" block. Outputs include pseudotime tables, gene–pseudotime correlations, and top-10 gene UMAP plots.

### 5. Differential expression

Run DESeq2-based differential expression analysis between two groups.

```r
source("5.DE_analysis/DE_analysis.R")
```

Set `paired = TRUE/FALSE` depending on experimental design. For paired analysis, ensure the metadata file contains the pairing column (default: `pair`). Edit `count_file`, `meta_file`, and `output_dir` in the "User settings" block.

### 6. Correlation

Compute Pearson and Spearman correlations of mean gene expression between two groups.

```r
source("6.correlation/correlation.R")
```

Edit `input_file`, `group_file`, and `output_dir` in the "User settings" block.


## License

Copyright (c) All authors of scBiopsy-seq. Released under the MIT License.
