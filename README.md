# DESeq2 Project

A differential gene expression analysis using [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) on the `airway` Bioconductor dataset (RNA-seq of airway smooth muscle cells treated with dexamethasone).

## Files

- `script_for_airwaypackage.R` — pulls the `airway` dataset from Bioconductor and writes out `counts_data.csv` and `sample_info.csv`.
- `counts_data.csv` — raw gene count matrix (genes x samples).
- `sample_info.csv` — sample metadata (treated / untreated, cell line).
- `DESeq2_project.R` — main analysis: builds the `DESeqDataSet`, filters low-count genes, runs `DESeq()`, and explores results (summary, contrasts, MA plot).
- `Rplot_DESeq2_project.pdf` — output MA plot from the analysis.
- `sample_qc_and_volcano.R` — sample QC (PCA + sample-distance heatmap) and a volcano plot annotated with gene symbols. Produces `pca_plot.png`, `sample_distance_heatmap.png`, `volcano_plot.png`, and `deseq2_results_annotated.csv`.

## Results

**Did the drug actually change the cells?** Yes. When samples are plotted by their overall gene activity (`pca_plot.png`), treated and untreated cells fall into two separate groups. So the drug's effect shows up clearly, just from the gene expression data.

**How many genes changed?** Out of about 18,400 genes we could test, **872 changed significantly** — 475 went up, 397 went down, after treatment (see `volcano_plot.png`).

**Do the results make sense?** Yes. Some of the top changed genes — PER1, DUSP1, ZBTB16, MAOA — are already known from other studies to respond to this exact type of drug (a steroid). Landing on the same genes independently is a good sign the analysis is actually working, not just producing random noise.

**One honest limitation:** part of the difference between samples comes from which person the cells were taken from, not only from the drug (visible in `sample_distance_heatmap.png`). A more advanced version of this analysis would account for that directly.

## Requirements

```r
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("DESeq2", "airway", "org.Hs.eg.db", "AnnotationDbi"))
install.packages(c("tidyverse", "pheatmap", "ggrepel", "RColorBrewer"))
```

## Usage

Run the scripts in order from this project's folder:

```r
source("script_for_airwaypackage.R")   # generates counts_data.csv / sample_info.csv
source("DESeq2_project.R")             # runs the DESeq2 analysis
source("sample_qc_and_volcano.R")      # PCA, sample-distance heatmap, volcano plot
```
