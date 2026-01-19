#heatmapChondro
# Load required libraries
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(dplyr)
library(readxl)

# 1) Read differential expression results and filter by padj
deg_results <- read_excel(
  "...differentialExpression_ThielChondr.xlsx"
) %>%
  filter(padj < 0.05)  # only significant DE genes

# 2) Define bona fide chondrogenic marker genes and subset
marker_genes <- c("SOX9", "COL2A1", "ACAN", "COMP", "COL10A1")
marker_deg_results <- deg_results %>%
  filter(gene_name %in% marker_genes) %>%
  distinct(ENSEMBL, gene_name, .keep_all = TRUE)

# 3) Load normalized (log) counts from sheet 3
tcounts <- read_excel(
  "...Counts_ThielRNA.xlsx",
  sheet = 3
) %>% as.data.frame()
stopifnot("ID" %in% names(tcounts))
rownames(tcounts) <- tcounts$ID; tcounts$ID <- NULL

# 4) Subset to chondrogenic + undifferentiated samples
cols_keep <- c(
  "WT_1","WT_2","WT_3",
  "KO_9","KO_75","KO_46",
  "CHONDR_WT_1","CHONDR_WT_2","CHONDR_WT_3",
  "CHONDR_KO_9","CHONDR_KO_75","CHONDR_KO_47"
)
log_counts <- tcounts[, cols_keep]

# 5) Filter rows to marker ENSEMBL IDs
gids <- intersect(marker_deg_results$ENSEMBL, rownames(log_counts))
if(length(gids) == 0) stop("No marker ENSEMBL IDs in count matrix after padj filtering")
marker_log_counts <- log_counts[gids, ]

# 6) Set rownames to gene symbols
genes <- marker_deg_results$gene_name[match(rownames(marker_log_counts), marker_deg_results$ENSEMBL)]
rownames(marker_log_counts) <- make.unique(genes)

# 7) Rename columns to concise sample names
colnames(marker_log_counts) <- c(
  "WT1","WT2","WT3",
  "KO1","KO2","KO3",
  "chWT1","chWT2","chWT3",
  "chKO1","chKO2","chKO3"
)

# 8) Annotation for samples
sample_annotation <- data.frame(
  Stage = rep(c("Undifferentiated","Chondrogenic"), each = 6),
  Genotype = c(
    rep("WT", 3), rep("KO", 3),  # Undifferentiated
    rep("WT", 3), rep("KO", 3)   # Chondrogenic
  ),
  row.names = colnames(marker_log_counts)
)

# 9) Define annotation colors
annotation_colors <- list(
  Stage    = c(Undifferentiated = "navy", Chondrogenic = "firebrick3"),
  Genotype = c(WT = "white", KO = "grey30")
)

# 10) Plot row-scaled heatmap of chondrogenic markers with larger fonts
pheatmap(
  marker_log_counts,
  scale             = "row",
  color             = colorRampPalette(c("navy","white","firebrick3"))(100),
  annotation_col    = sample_annotation,
  annotation_colors = annotation_colors,
  cluster_rows      = TRUE,
  cluster_cols      = FALSE,
  show_rownames     = TRUE,
  show_colnames     = TRUE,
  fontsize          = 16,
  fontsize_row      = 16,
  fontsize_col      = 16,
  cellwidth         = 25,
  cellheight        = 25,
  border_color      = "black",
  main              = "Chondrogenic Marker Expression"
)

#heatmapOsteo
# 1) Read DE results for osteogenesis and filter by padj
deg_results <- read_excel(
  "...differentialExpression_ThielOsteo.xlsx"
) %>%
  filter(padj < 0.05)

# 2) Define osteogenic markers and subset
marker_genes <- c("RUNX2", "SP7", "ALPL", "BGLAP", "SPP1")
marker_deg_results <- deg_results %>%
  filter(gene_name %in% marker_genes) %>%
  distinct(ENSEMBL, gene_name, .keep_all = TRUE)

# 3) Load counts (already in tcounts)
# 4) Subset to undifferentiated and osteogenic samples
cols_keep <- c(
  "WT_1","WT_2","WT_3",
  "KO_9","KO_75","KO_46",
  "OsteoLAMA5_WT_1","OsteoLAMA5_WT_2","OsteoLAMA5_WT_3",
  "OsteoLAMA5_KO_9","OsteoLAMA5_KO_75","OsteoLAMA5_KO_46"
)
log_counts <- tcounts[, cols_keep]

# 5) Filter rows to marker IDs
gids <- intersect(marker_deg_results$ENSEMBL, rownames(log_counts))
if(length(gids) == 0) stop("No osteogenic marker ENSEMBL IDs in count matrix after padj filtering")
marker_log_counts <- log_counts[gids, ]

# 6) Rename rows to gene symbols
gene_symbols <- marker_deg_results$gene_name[match(rownames(marker_log_counts), marker_deg_results$ENSEMBL)]
rownames(marker_log_counts) <- make.unique(gene_symbols)

# 7) Rename columns concisely
colnames(marker_log_counts) <- c(
  "WT1","WT2","WT3",
  "KO1","KO2","KO3",
  "oWT1","oWT2","oWT3",
  "oKO1","oKO2","oKO3"
)

# 8) Sample annotation
sample_annotation <- data.frame(
  Stage    = rep(c("Undifferentiated","Osteogenic"), each = 6),
  Genotype = c(
    rep("WT", 3), rep("KO", 3),
    rep("WT", 3), rep("KO", 3)
  ),
  row.names = colnames(marker_log_counts)
)
annotation_colors <- list(
  Stage    = c(Undifferentiated = "navy", Osteogenic = "goldenrod3"),
  Genotype = c(WT = "white", KO = "grey30")
)

# 10) Plot heatmap of osteogenic markers with larger fonts
pheatmap(
  marker_log_counts,
  scale             = "row",
  color             = colorRampPalette(c("navy","white","goldenrod3"))(100),
  annotation_col    = sample_annotation,
  annotation_colors = annotation_colors,
  cluster_rows      = TRUE,
  cluster_cols      = FALSE,
  show_rownames     = TRUE,
  show_colnames     = TRUE,
  fontsize          = 16,
  fontsize_row      = 16,
  fontsize_col      = 16,
  cellwidth         = 25,
  cellheight        = 25,
  border_color      = "black",
  main              = "Osteogenic Marker Expression"
)
