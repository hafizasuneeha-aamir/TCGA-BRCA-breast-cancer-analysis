# ============================================================
# Differential Gene Expression Analysis of Breast Cancer
# TCGA-BRCA: Tumour vs Normal Tissue
# Author: Hafiza Suneeha Aamir
# ============================================================

# ---- 1. Install and load required packages ----
if (!require('BiocManager', quietly = TRUE))
  install.packages('BiocManager')

BiocManager::install(c('TCGAbiolinks', 'DESeq2', 'EnhancedVolcano',
                        'clusterProfiler', 'org.Hs.eg.db', 'ComplexHeatmap'))
install.packages(c('ggplot2', 'dplyr', 'pheatmap'))

library(TCGAbiolinks)
library(DESeq2)
library(dplyr)
library(EnhancedVolcano)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)

# ---- 2. Download TCGA-BRCA RNA-seq data ----
options(timeout = 600)

query <- GDCquery(
  project = 'TCGA-BRCA',
  data.category = 'Transcriptome Profiling',
  data.type = 'Gene Expression Quantification',
  workflow.type = 'STAR - Counts'
)

GDCdownload(query, files.per.chunk = 50)
data <- GDCprepare(query)

dim(data)  # 60660 genes x 1231 samples

# ---- 3. Label samples as tumour or normal ----
counts  <- assay(data, "unstranded")
coldata <- as.data.frame(colData(data))

table(coldata$sample_type)

coldata$condition <- case_when(
  coldata$sample_type == "Primary Tumor"       ~ "tumour",
  coldata$sample_type == "Solid Tissue Normal"  ~ "normal",
  TRUE ~ NA_character_
)

keep    <- !is.na(coldata$condition)
counts  <- counts[, keep]
coldata <- coldata[keep, ]

coldata$condition <- factor(coldata$condition, levels = c("normal", "tumour"))
table(coldata$condition)  # normal: 113, tumour: 1111

# ---- 4. Differential expression with DESeq2 ----
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = coldata,
  design    = ~ condition
)

dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)

res <- results(dds, contrast = c('condition', 'tumour', 'normal'))
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)

head(res_df[order(res_df$padj), ])

write.csv(res_df, "results/results_table.csv", row.names = FALSE)

# ---- 5. Volcano plot ----
EnhancedVolcano(res_df,
  lab = res_df$gene,
  x = 'log2FoldChange',
  y = 'padj',
  pCutoff = 0.05,
  FCcutoff = 1,
  title = 'Tumour vs Normal - Differential Gene Expression',
  subtitle = 'TCGA-BRCA',
  legendPosition = 'right'
)
# Export manually or with ggsave() as figures/volcano_plot.png

# ---- 6. Heatmap of top 50 DEGs ----
sig_genes <- res_df %>%
  filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 1) %>%
  arrange(padj) %>%
  head(50)

vsd <- vst(dds, blind = FALSE)
mat <- assay(vsd)[sig_genes$gene, ]
mat_scaled <- t(scale(t(mat)))

annotation <- data.frame(
  condition = coldata$condition,
  row.names = colnames(mat)
)

pheatmap(mat_scaled,
  annotation_col = annotation,
  show_rownames = TRUE,
  show_colnames = FALSE,
  main = 'Top 50 DEGs: Tumour vs Normal'
)
# Export manually as figures/heatmap.png

# ---- 7. GO enrichment analysis (biological process) ----
up_genes <- res_df %>%
  filter(!is.na(padj), padj < 0.05, log2FoldChange > 1) %>%
  pull(gene)

up_genes_clean <- sub("\\..*", "", up_genes)

gene_ids <- bitr(up_genes_clean,
  fromType = 'ENSEMBL',
  toType = 'ENTREZID',
  OrgDb = org.Hs.eg.db)

ego <- enrichGO(
  gene = gene_ids$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = 'BP',
  pAdjustMethod = 'BH',
  pvalueCutoff = 0.05,
  readable = TRUE
)

dotplot(ego, showCategory = 20, title = 'GO Biological Process Enrichment')
# Export manually as figures/go_enrichment.png

# ---- 8. Gene name lookup for top hits ----
top_genes <- res_df %>%
  filter(!is.na(padj)) %>%
  arrange(padj) %>%
  head(10) %>%
  mutate(ensembl_clean = sub("\\..*", "", gene))

top_symbols <- bitr(top_genes$ensembl_clean,
  fromType = 'ENSEMBL',
  toType = 'SYMBOL',
  OrgDb = org.Hs.eg.db)

top_symbols

# ============================================================
# Summary of results:
#   12,304 significantly differentially expressed genes
#   (8,186 upregulated, 4,118 downregulated)
#   Top hits: MMP11, COL10A1, NEK2, COL11A1, KIF4A (up)
#             MYOM1, GPAM, FHL1, TNS1 (down)
#   GO enrichment: cell division / chromosome segregation dominant
# ============================================================