library(ggplot2)
library(dplyr)
library(ggrepel)
library(readxl)

# ------------------- General Volcano Plot Function ------------------- #

create_volcano_plot <- function(file_path, condition_label, marker_genes, colors, title_text = NULL) {
  
  # Load and clean table
  deg_table <- read_xlsx(file_path)
  deg_table_clean <- deg_table %>%
    filter(!is.na(log2FoldChange), !is.na(padj))
  deg_table_clean$log2FoldChange <- -deg_table_clean$log2FoldChange
  
  # Thresholds
  pval_threshold <- 0.05
  logfc_threshold <- 1
  
  # Categorize by significance and direction
  deg_table_clean <- deg_table_clean %>%
    mutate(
      Significance = case_when(
        padj < pval_threshold & abs(log2FoldChange) > logfc_threshold ~ "Significant",
        TRUE ~ "Not significant"
      ),
      Direction = case_when(
        padj < pval_threshold & log2FoldChange > logfc_threshold ~ paste0(condition_label, " up"),
        padj < pval_threshold & log2FoldChange < -logfc_threshold ~ paste0(condition_label, " down"),
        TRUE ~ "Not significant"
      ),
      Color = Direction,
      Marker = ifelse(gene_name %in% marker_genes, gene_name, NA)
    )
  
  # Filter extreme outliers
  filtered_deg_table <- deg_table_clean %>%
    filter(between(log2FoldChange, -10, 10))
  
  # Plot
  volcano_plot <- ggplot(filtered_deg_table, aes(x = log2FoldChange, y = -log10(padj), color = Color)) +
    geom_point(alpha = 0.85, size = 3) +
    scale_color_manual(values = colors) +
    geom_vline(xintercept = c(-logfc_threshold, logfc_threshold), 
               linetype = "dotted", color = "grey40", linewidth = 0.6) +
    geom_hline(yintercept = -log10(pval_threshold), 
               linetype = "dotted", color = "grey40", linewidth = 0.6) +
    geom_point(data = filtered_deg_table %>% filter(!is.na(Marker)), 
               size = 4.2, color = "black", shape = 21, fill = "white", stroke = 1.2) +
    geom_label_repel(
      data = filtered_deg_table %>%
        filter(!is.na(Marker) & padj < pval_threshold & abs(log2FoldChange) > logfc_threshold),
      aes(label = Marker),
      color = "black", 
      size = 6,
      fill = "white",
      fontface = "bold.italic",
      box.padding = 0.3, 
      point.padding = 0.3,
      max.overlaps = 20
    ) +
    labs(
      title = title_text,
      x = expression(log[2]~"Fold Change"),
      y = expression(-log[10]~"(p-adj.)"),
      color = NULL
    ) +
    theme_classic(base_size = 16) +
    theme(
      legend.position = "none",
      plot.title      = element_text(size = 18, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 16, face = "bold"),
      axis.text = element_text(size = 14),
      axis.line = element_line(size = 0.8),
      panel.grid = element_blank()
    )
  
  return(volcano_plot)
}

# ------------------- Plot CDH1 ------------------- #

chondro_path <- "differentialExpression_Thiel_Chond.xlsx..."
chondro_markers <- c(
  "WNT1","WNT2","WNT2B","WNT3","WNT3A","WNT4","WNT5A","WNT5B",
  "WNT6","WNT7A","WNT7B","WNT8A","WNT8B","WNT9A","WNT9B",
  "WNT10A","WNT10B","WNT11","WNT16",
  "FZD1","FZD2","FZD3","FZD4","FZD5","FZD6","FZD7","FZD8","FZD9","FZD10",
  "LRP5","LRP6","ROR1","ROR2","RYK",
  "CTNNB1","APC","AXIN1","AXIN2","GSK3B",
  "CSNK1A1","CSNK1D","CSNK1E",
  "DVL1","DVL2","DVL3",
  "TCF7","TCF7L1","TCF7L2","LEF1",
  "BCL9","BCL9L","PYGO1","PYGO2",
  "DKK1","DKK2","DKK3","DKK4",
  "SFRP1","SFRP2","SFRP3","SFRP4","SFRP5",
  "WIF1","NKD1","NKD2",
  "MYC","CCND1","LGR5","JUN","FOSL1","MMP7"
)

chondro_colors <- c("LAMA5 KO up" = "#FF6666", 
                    "LAMA5 KO down" = "#66A3FF", 
                    "Not significant" = "grey80")

volcano_chondro <- create_volcano_plot(
  file_path = chondro_path,
  condition_label = "LAMA5 KO",
  marker_genes = chondro_markers,
  colors = chondro_colors,
  title_text = "LAMA5 KO vs. WT"
)

print(volcano_chondro)




library(dplyr)
library(readxl)
library(pheatmap)
library(RColorBrewer)

# 1) Load DEG results (filtered for significance)
deg_file <- "differentialExpression_Thiel_Chond.xlsx..."
deg_results <- read_excel(deg_file) %>% filter(padj < 0.05)

# 2) Define your WNT gene set (same as before)
wnt_genes <- c(
  "WNT1","WNT2","WNT2B","WNT3","WNT3A","WNT4","WNT5A","WNT5B",
  "WNT6","WNT7A","WNT7B","WNT8A","WNT8B","WNT9A","WNT9B",
  "WNT10A","WNT10B","WNT11","WNT16",
  "FZD1","FZD2","FZD3","FZD4","FZD5","FZD6","FZD7","FZD8","FZD9","FZD10",
  "LRP5","LRP6","ROR1","ROR2","RYK",
  "CTNNB1","APC","AXIN1","AXIN2","GSK3B",
  "CSNK1A1","CSNK1D","CSNK1E",
  "DVL1","DVL2","DVL3",
  "TCF7","TCF7L1","TCF7L2","LEF1",
  "BCL9","BCL9L","PYGO1","PYGO2",
  "DKK1","DKK2","DKK3","DKK4",
  "SFRP1","SFRP2","SFRP3","SFRP4","SFRP5",
  "WIF1","NKD1","NKD2",
  "MYC","CCND1","LGR5","JUN","FOSL1","MMP7"
)

# 3) Subset significant DEGs to WNT genes only
wnt_deg <- deg_results %>% filter(gene_name %in% wnt_genes)

# 4) Load normalized (log) counts
counts_file <- "Counts_ThielRNA.xlsx..."
log_counts <- read_excel(counts_file, sheet = 3) %>% as.data.frame()

stopifnot("ID" %in% names(log_counts))
rownames(log_counts) <- log_counts$ID
log_counts$ID <- NULL

# 5) Subset to WNT genes present in counts and filtered DEGs
gids <- intersect(wnt_deg$ENSEMBL, rownames(log_counts))
if(length(gids) == 0) stop("No WNT ENSEMBL IDs found in counts matrix after filtering")

wnt_log_counts <- log_counts[gids, ]

# 6) Replace rownames with gene symbols for clarity
gene_symbols <- wnt_deg$gene_name[match(rownames(wnt_log_counts), wnt_deg$ENSEMBL)]
rownames(wnt_log_counts) <- make.unique(gene_symbols)

# 7) Select relevant samples (edit to your experiment design)
samples_keep <- c(
  "WT_1","WT_2","WT_3",
  "KO_9","KO_75","KO_46"
  # Add more if you want (e.g. chondrogenic or osteogenic samples)
)
wnt_log_counts <- wnt_log_counts[, samples_keep]

# 8) Define sample annotations
sample_annotation <- data.frame(
  Genotype = c(rep("WT", 3), rep("KO", 3)),
  row.names = samples_keep
)

annotation_colors <- list(
  Genotype = c(WT = "white", KO = "grey30")
)

# 9) Plot heatmap
pheatmap(
  wnt_log_counts,
  scale             = "row",
  color             = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  annotation_col    = sample_annotation,
  annotation_colors = annotation_colors,
  cluster_rows      = TRUE,
  cluster_cols      = TRUE,
  show_rownames     = TRUE,
  show_colnames     = TRUE,
  fontsize          = 14,
  fontsize_row      = 14,
  fontsize_col      = 14,
  cellwidth         = 25,
  cellheight        = 25,
  border_color      = "black",
  main              = "Expression Heatmap of Significant WNT Genes"
)

# Step 7: Rename columns concisely (samples)
colnames(wnt_log_counts) <- c("WT1", "WT2", "WT3", "KO1", "KO2", "KO3")  # adjust as needed

# Step 7.5: Transpose matrix (samples as rows, genes as columns)
wnt_log_counts_t <- t(wnt_log_counts)

# Step 8: Define sample annotations for the transposed matrix
# Now rownames = samples
sample_annotation <- data.frame(
  Genotype = c(rep("WT", 3), rep("KO", 3)),
  row.names = rownames(wnt_log_counts_t)
)

annotation_colors <- list(
  Genotype = c(WT = "white", KO = "grey30")
)

# Step 9: Plot heatmap with flipped orientation
pheatmap(
  wnt_log_counts_t,
  scale             = "column",   # scale by gene now (columns)
  color             = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  annotation_row    = sample_annotation,   # samples are rows
  annotation_colors = annotation_colors,
  cluster_rows      = TRUE,
  cluster_cols      = TRUE,
  show_rownames     = TRUE,
  show_colnames     = TRUE,
  fontsize          = 14,
  fontsize_row      = 14,
  fontsize_col      = 14,
  cellwidth         = 25,
  cellheight        = 25,
  border_color      = "black",
  main              = "Expression Heatmap of Significant WNT Genes"
)