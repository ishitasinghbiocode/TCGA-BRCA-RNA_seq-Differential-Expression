# ============================================================
# TCGA-BRCA Differential Gene Expression Analysis
# Tumour vs Normal Breast Tissue
# ============================================================

# -------------------------
# 1. Load required packages
# -------------------------

library(TCGAbiolinks)
library(SummarizedExperiment)
library(DESeq2)
library(EnhancedVolcano)
library(ggplot2)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

# -------------------------
# 2. Define TCGA project
# -------------------------

project <- "TCGA-BRCA"

# -------------------------
# 3. Query TCGA-BRCA data
# -------------------------

query <- GDCquery(
  project = project,
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor", "Solid Tissue Normal")
)

# -------------------------
# 4. Select balanced samples
# 50 tumour + 50 normal
# -------------------------

results <- getResults(query)

set.seed(123)

tumor_barcodes <- sample(
  results$cases[results$sample_type == "Primary Tumor"],
  50
)

normal_barcodes <- sample(
  results$cases[results$sample_type == "Solid Tissue Normal"],
  50
)

selected_barcodes <- c(
  tumor_barcodes,
  normal_barcodes
)

# -------------------------
# 5. Query selected samples
# -------------------------

query_100 <- GDCquery(
  project = project,
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  barcode = selected_barcodes
)

# -------------------------
# 6. Download and prepare data
# -------------------------

GDCdownload(query_100)

data <- GDCprepare(query_100)

# -------------------------
# 7. Extract count matrix
# -------------------------

counts <- assay(data, "unstranded")

coldata <- as.data.frame(colData(data))

condition <- factor(
  ifelse(
    coldata$sample_type == "Primary Tumor",
    "tumour",
    "normal"
  ),
  levels = c("normal", "tumour")
)

# -------------------------
# 8. Create DESeq2 object
# -------------------------

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = data.frame(condition = condition),
  design = ~condition
)

# -------------------------
# 9. Filter low-count genes
# -------------------------

keep <- rowSums(counts(dds) >= 10) >= 10

dds <- dds[keep, ]

# -------------------------
# 10. Differential expression
# -------------------------

dds <- DESeq(dds)

res <- results(
  dds,
  contrast = c(
    "condition",
    "tumour",
    "normal"
  )
)

res_df <- as.data.frame(res)

res_df$gene_id <- rownames(res_df)

res_df <- res_df[
  order(res_df$padj),
]

# -------------------------
# 11. Identify significant DEGs
# padj < 0.05 and |log2FC| > 1
# -------------------------

sig <- subset(
  res_df,
  padj < 0.05 &
    abs(log2FoldChange) > 1
)

# -------------------------
# 12. Create output folders
# -------------------------

dir.create(
  "results",
  showWarnings = FALSE
)

dir.create(
  "figures",
  showWarnings = FALSE
)

# -------------------------
# 13. Save DEG results
# -------------------------

write.csv(
  res_df,
  "results/DEG_results.csv",
  row.names = FALSE
)

write.csv(
  sig,
  "results/significant_DEGs.csv",
  row.names = FALSE
)

# -------------------------
# 14. Volcano plot
# -------------------------

volcano_df <- res_df[
  !is.na(res_df$padj),
]

top_labels <- head(
  res_df[order(res_df$padj), ],
  20
)

volcano_plot <- EnhancedVolcano(
  volcano_df,
  lab = ifelse(
    volcano_df$gene_id %in%
      top_labels$gene_id,
    volcano_df$gene_id,
    ""
  ),
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.05,
  FCcutoff = 1,
  title = "TCGA-BRCA: Tumour vs Normal",
  subtitle = "Differential Gene Expression",
  caption = "Significant: adjusted p-value < 0.05 and |log2FC| > 1"
)

ggsave(
  "figures/volcano_plot.png",
  plot = volcano_plot,
  width = 10,
  height = 7,
  dpi = 300
)

# -------------------------
# 15. Variance stabilizing transformation
# -------------------------

vsd <- vst(
  dds,
  blind = FALSE
)

# -------------------------
# 16. Heatmap of top 50 DEGs
# -------------------------

top_genes <- head(
  sig$gene_id,
  50
)

mat <- assay(vsd)[
  top_genes,
  ,
  drop = FALSE
]

mat <- t(
  scale(
    t(mat)
  )
)

annotation_col <- data.frame(
  condition = condition
)

rownames(annotation_col) <-
  colnames(mat)

pheatmap(
  mat,
  annotation_col = annotation_col,
  show_rownames = FALSE,
  main = "Top 50 Differentially Expressed Genes",
  filename = "figures/heatmap_top50.png",
  width = 10,
  height = 9
)

# -------------------------
# 17. Prepare genes for GO analysis
# -------------------------

sig_ids <- sub(
  "\\..*$",
  "",
  sig$gene_id
)

gene_map <- bitr(
  sig_ids,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

# -------------------------
# 18. GO Biological Process
# enrichment
# -------------------------

ego <- enrichGO(
  gene = unique(
    gene_map$ENTREZID
  ),
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)

# -------------------------
# 19. GO dot plot
# -------------------------

p_go <- dotplot(
  ego,
  showCategory = 15,
  title = "GO Biological Process Enrichment"
)

ggsave(
  "figures/GO_BP_dotplot.png",
  plot = p_go,
  width = 10,
  height = 8,
  dpi = 300
)

# -------------------------
# 20. Save GO results
# -------------------------

go_df <- as.data.frame(ego)

write.csv(
  go_df,
  "results/GO_BP_enrichment.csv",
  row.names = FALSE
)

# -------------------------
# 21. Save DESeq2 object
# -------------------------

saveRDS(
  dds,
  "results/BRCA_DESeq2.rds"
)

# -------------------------
# 22. Save session information
# -------------------------

writeLines(
  capture.output(
    sessionInfo()
  ),
  "results/sessionInfo.txt"
)

# ============================================================
# END OF ANALYSIS
# ============================================================