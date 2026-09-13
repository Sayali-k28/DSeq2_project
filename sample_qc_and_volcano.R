# Sample QC (PCA + sample-distance heatmap) and a volcano plot with real
# gene symbols, for the dexamethasone-treated vs. untreated airway smooth
# muscle cell DESeq2 analysis.
#
# This rebuilds the same DESeqDataSet as DESeq2_project.R (kept self-contained
# so it runs on its own) and adds the diagnostic + result plots that script
# doesn't produce.

library(DESeq2)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(ggrepel)
library(AnnotationDbi)
library(org.Hs.eg.db)

#step 0: rebuild the same DESeqDataSet as DESeq2_project.R -----------
counts_data <- read.csv('counts_data.csv')
colData <- read.csv('sample_info.csv')

dds <- DESeqDataSetFromMatrix(countData = counts_data,
                               colData = colData,
                               design = ~ dexamethasone)

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]
dds$dexamethasone <- relevel(dds$dexamethasone, ref = "untreated")
dds <- DESeq(dds)
res <- results(dds)

#step 1a: PCA plot -----------
# variance-stabilizing transform: makes counts roughly homoskedastic so
# Euclidean-distance-based methods like PCA behave well (raw counts have
# variance that scales with the mean, which would distort this).
vsd <- vst(dds, blind = TRUE)

pca_data <- plotPCA(vsd, intgroup = c("dexamethasone", "cellLine"), returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))

p_pca <- ggplot(pca_data, aes(PC1, PC2, color = dexamethasone, shape = cellLine)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", percent_var[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var[2], "% variance")) +
  ggtitle("PCA of samples after variance-stabilizing transform") +
  theme_bw(base_size = 13)

ggsave("pca_plot.png", p_pca, width = 6.5, height = 5, dpi = 150)

#step 1b: sample-to-sample distance heatmap -----------
sample_dists <- dist(t(assay(vsd)))
dist_matrix <- as.matrix(sample_dists)
rownames(dist_matrix) <- paste(vsd$dexamethasone, vsd$cellLine, sep = " - ")
colnames(dist_matrix) <- NULL

png("sample_distance_heatmap.png", width = 1400, height = 1200, res = 200)
pheatmap(dist_matrix,
         clustering_distance_rows = sample_dists,
         clustering_distance_cols = sample_dists,
         col = colorRampPalette(rev(brewer.pal(9, "Blues")))(255),
         main = "Sample-to-sample distances (VST-transformed counts)")
dev.off()

cat("Wrote pca_plot.png and sample_distance_heatmap.png\n")

#step 3: volcano plot with gene symbols -----------
res_df <- as.data.frame(res)
res_df$ensembl_id <- rownames(res_df)
res_df <- res_df[!is.na(res_df$padj), ]

# map Ensembl -> gene symbol; not every gene has one, that's expected
res_df$symbol <- mapIds(org.Hs.eg.db,
                         keys = res_df$ensembl_id,
                         column = "SYMBOL",
                         keytype = "ENSEMBL",
                         multiVals = "first")

res_df$significance <- case_when(
  res_df$padj < 0.05 & res_df$log2FoldChange > 1  ~ "Up in treated",
  res_df$padj < 0.05 & res_df$log2FoldChange < -1 ~ "Down in treated",
  TRUE ~ "Not significant"
)

# label the 10 most significant genes that also have a known symbol
top_labels <- res_df %>%
  filter(!is.na(symbol), significance != "Not significant") %>%
  arrange(padj) %>%
  head(10)

p_volcano <- ggplot(res_df, aes(log2FoldChange, -log10(padj), color = significance)) +
  geom_point(alpha = 0.5, size = 1.2) +
  scale_color_manual(values = c("Up in treated" = "#B2182B",
                                 "Down in treated" = "#2166AC",
                                 "Not significant" = "grey70")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  geom_text_repel(data = top_labels, aes(label = symbol),
                   color = "black", size = 3.3, max.overlaps = 20) +
  labs(title = "Dexamethasone-treated vs. untreated airway smooth muscle cells",
       subtitle = "Thresholds: padj < 0.05, |log2FC| > 1",
       x = "log2 fold change", y = "-log10(adjusted p-value)", color = NULL) +
  theme_bw(base_size = 13)

ggsave("volcano_plot.png", p_volcano, width = 7.5, height = 6, dpi = 150)

cat("Wrote volcano_plot.png\n")
cat("\nTop 10 most significant genes with a known symbol:\n")
print(top_labels[, c("symbol", "log2FoldChange", "padj")])

# save the annotated table too, useful for the README + future steps
write.csv(res_df, "deseq2_results_annotated.csv", row.names = FALSE)
cat("\nWrote deseq2_results_annotated.csv\n")
