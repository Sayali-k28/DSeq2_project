# DESeq2 Project

A differential gene expression analysis using [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) on the `airway` Bioconductor dataset (RNA-seq of airway smooth muscle cells treated with dexamethasone).

## Files

- `script_for_airwaypackage.R` — pulls the `airway` dataset from Bioconductor and writes out `counts_data.csv` and `sample_info.csv`.
- `counts_data.csv` — raw gene count matrix (genes x samples).
- `sample_info.csv` — sample metadata (treated / untreated, cell line).
- `DESeq2_project.R` — main analysis: builds the `DESeqDataSet`, filters low-count genes, runs `DESeq()`, and explores results (summary, contrasts, MA plot).
- `Rplot_DESeq2_project.pdf` — output MA plot from the analysis.

## Requirements

```r
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("DESeq2", "airway"))
install.packages("tidyverse")
```

## Usage

Run the scripts in order from this project's folder:

```r
source("script_for_airwaypackage.R")  # generates counts_data.csv / sample_info.csv
source("DESeq2_project.R")            # runs the DESeq2 analysis
```
