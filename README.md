# Differential Gene Expression Analysis of Breast Cancer (TCGA-BRCA)

Analysis of RNA-seq data from The Cancer Genome Atlas (TCGA-BRCA) comparing
gene expression in breast tumor tissue versus matched normal tissue, using
DESeq2 for differential expression and GO enrichment for pathway analysis.

## Summary

- **Dataset:** TCGA-BRCA, 1,231 samples (1,111 tumor, 113 normal)
- **Method:** DESeq2 differential expression, clusterProfiler GO enrichment
- **Key finding:** 12,304 significantly differentially expressed genes
  (8,186 upregulated, 4,118 downregulated). Top hits include MMP11 and
  COL11A1, both independently validated breast cancer biomarkers in
  published literature. GO enrichment of upregulated genes showed strong
  enrichment for cell division and chromosome segregation processes.

## Repository structure

- `code/` — R script used for the full analysis
- `figures/` — volcano plot, heatmap, GO enrichment plot
- `results/` — full differential expression results table (CSV)
- `report/` — full written report (background, methods, results, discussion)

## Tools used

R, DESeq2, TCGAbiolinks, clusterProfiler, EnhancedVolcano, pheatmap

## Author

Hafiza Suneeha Aamir
