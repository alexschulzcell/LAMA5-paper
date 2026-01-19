library(readxl)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(biomaRt)
library(ReactomePA)
library(RColorBrewer)


deg_table <- read_xlsx("differentialExpression_Thiel_Chond.xlsx...")
colnames(deg_table) <- c("ENSEMBL", "gene_name","gene_type", "gene_description", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "hyperlink")
deg_table_clean <- deg_table %>%
  filter(!is.na(log2FoldChange) & !is.na(padj))
pval_threshold <- 0.05
logfc_threshold <- 1
deg_table_clean <- deg_table_clean %>%
  filter(padj < pval_threshold & abs(log2FoldChange) > logfc_threshold)
# Set up biomaRt
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Map Gene Symbols to Entrez Gene IDs
entrez_mapping <- getBM(
  attributes = c("external_gene_name", "entrezgene_id"),
  filters = "external_gene_name",
  values = deg_table_clean$gene_name,
  mart = ensembl
)
# Perform GO enrichment analysis for Biological Process (BP)
go_enrichment_BP <- enrichGO(
  gene = entrez_mapping$entrezgene_id,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",  # Biological Process
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = TRUE
)

# Perform GO enrichment analysis for Cellular Component (CC)
go_enrichment_CC <- enrichGO(
  gene = entrez_mapping$entrezgene_id,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "CC",  # Cellular Component
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = TRUE
)

# Perform GO enrichment analysis for Molecular Function (MF)
go_enrichment_MF <- enrichGO(
  gene = entrez_mapping$entrezgene_id,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "MF",  # Molecular Function
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = TRUE
)
# Reactome Pathway Enrichment Analysis
reactome_enrichment <- enrichPathway(
  gene = entrez_mapping$entrezgene_id,
  organism = 'human',
  pvalueCutoff = 0.05
)
# GO Enrichment Plot for Biological Process (BP)
dotplot(go_enrichment_BP, showCategory = 10) +
  ggtitle("Gene Ontology Enrichment - Biological Processes") +
  theme_minimal(base_size = 14) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
        plot.title = element_text(size = 16, face = "bold"),  # Title styling
        axis.title = element_text(size = 14))

# GO Enrichment Plot for Cellular Components (CC)
dotplot(go_enrichment_CC, showCategory = 10) + 
  ggtitle("Gene Ontology Enrichment - Cellular Components") + 
  theme_minimal(base_size = 14) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 14))

# GO Enrichment Plot for Molecular Functions (MF)
dotplot(go_enrichment_MF, showCategory = 10) + 
  ggtitle("Gene Ontology Enrichment - Molecular Functions") + 
  theme_minimal(base_size = 14) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 14))
library(RColorBrewer)
library(enrichplot)
library(ggplot2)
library(dplyr)

# --- prepare the named foldChange vector correctly ---
fc_df <- deg_table_clean %>%
  inner_join(entrez_mapping, 
             by = c("gene_name" = "external_gene_name")) %>%
  dplyr::select(entrezgene_id, log2FoldChange)

fc <- with(fc_df, setNames(log2FoldChange, as.character(entrezgene_id)))
# invert the sign
fc <- -fc

head(fc[order(fc)])  # Should show strongly downregulated genes (negative)
head(fc[order(-fc)]) # Should show strongly upregulated genes (positive)

p <- cnetplot(
  reactome_enrichment,
  showCategory       = 8,
  foldChange         = fc,
  colorEdge          = FALSE,
  circular           = FALSE,
  node_label         = "category",
  edge_width         = 0.1,
  layout             = "kk",
  shadowtext         = "category"
)

# Manually override label size in layer 4 (category labels)
p$layers[[4]]$aes_params$size <- 5  # or whatever you prefer, like 8 or 10

# Add your scales and themes as before
p +
  scale_color_gradient2(
    low      = "navy", 
    mid      = "white",
    high     = "firebrick3", 
    midpoint = 0,
    limits   = range(fc, na.rm = TRUE),
    name     = "log2FC"
  ) +
  theme_void() +
  theme(
    legend.position = "right",
    legend.title    = element_text(face = "bold", size = 14),
    legend.text     = element_text(size = 12),
    plot.title      = element_text(size = 18, face = "bold", hjust = 0.5)
  ) +
  labs(title = "Reactome Pathway Enrichment")

