# --- Load libraries ---

library(readxl)
library(WGCNA)
library(clusterProfiler)
library(org.Hs.eg.db)
library(biomaRt)
library(DOSE)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(dplyr)
library(enrichplot)

# Enable multithreading for WGCNA
allowWGCNAThreads()

# --- Load counts data ---

# Assumes: first column = gene IDs (Ensembl), remaining columns = counts per sample.
counts_data <- read_excel("...Counts_ThielRNA.xlsx", sheet = 3)

gene_ids <- counts_data[[1]]         # Gene IDs (Ensembl)
counts_data_counts_only <- counts_data[, -1]  # Counts only

# --- Prepare expression matrix for WGCNA ---

# Your Excel file likely has genes in rows and samples in columns.
# For WGCNA, the standard input is a data matrix with *samples as rows* and *genes as columns*.
# Therefore, we transpose the counts matrix and assign gene IDs to the columns.
counts_matrix <- as.data.frame(t(counts_data_counts_only))
colnames(counts_matrix) <- gene_ids          # Now, columns = gene IDs
# (Row names are currently sample names from the Excel header)

cat("Samples:", nrow(counts_matrix), "Genes:", ncol(counts_matrix), "/n")

# --- Filter lowly expressed genes ---
# Now, genes are in columns so use colSums().
keep_genes <- colSums(counts_matrix > 1) >= 4
counts_matrix_filtered <- counts_matrix[, keep_genes]

# --- Log2 transformation ---
counts_matrix_normalized <- log2(counts_matrix_filtered + 1)

# --- Rename sample names (rownames) ---

samples <- rownames(counts_matrix_normalized)
samples <- gsub("OsteoLAMA5_WT_", "oWT", samples)
samples <- gsub("OsteoLAMA5_KO_", "oKO", samples)
samples <- gsub("CHONDR_WT_",   "chWT", samples)
samples <- gsub("CHONDR_KO_",   "chKO", samples)
samples <- gsub("KO9",          "KO1", samples)
samples <- gsub("KO75",         "KO2", samples)
samples <- gsub("KO46",         "KO3", samples)
samples <- gsub("KO_9",         "KO1", samples)
samples <- gsub("KO_75",        "KO2", samples)
samples <- gsub("KO_46",        "KO3", samples)
samples <- gsub("WT_1",         "WT1", samples)
samples <- gsub("WT_2",         "WT2", samples)
samples <- gsub("WT_3",         "WT3", samples)
rownames(counts_matrix_normalized) <- samples

# --- Pick soft-thresholding power ---
powers <- c(1:20)
sft <- pickSoftThreshold(counts_matrix_normalized, powerVector = powers, verbose = 5)

# Plot scale-free topology model fit
plot(sft$fitIndices[, 1], sft$fitIndices[, 2], type = "n",
     xlab = "Soft Threshold (Power)",
     ylab = "Scale-free Topology Model Fit (R^2)",
     main = "Scale Independence")
text(sft$fitIndices[, 1], sft$fitIndices[, 2], labels = powers, col = "red")
abline(h = 0.8, col = "blue", lty = 2)  # Threshold line

# Choose a power based on the plot (adjust as needed)
chosen_power <- 9

cor <- WGCNA::cor
# --- Build network and detect modules ---

# WGCNA expects data with samples as rows and genes as columns.
net <- blockwiseModules(
  counts_matrix_normalized,
  power = chosen_power,
  TOMType = "unsigned",
  minModuleSize = 30,
  reassignThreshold = 0,
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  verbose = 3
)

# Convert numeric labels to colors.
module_colors <- labels2colors(net$colors)
# Assign gene IDs (which are in the column names) to the module_colors vector.
names(module_colors) <- colnames(counts_matrix_normalized)

# --- Trait Analysis ---

# Sample-level traits. Here, use the sample names (rownames) of counts_matrix_normalized.
sample_status <- ifelse(grepl("KO", rownames(counts_matrix_normalized)), "KO", "WT")
traits <- data.frame(KO_WT_status = factor(sample_status, levels = c("WT", "KO")))
stopifnot(nrow(traits) == nrow(counts_matrix_normalized))

# Recalculate module eigengenes (MEs). WGCNA calculates MEs from the expression data.
MEs <- moduleEigengenes(counts_matrix_normalized, colors = module_colors)$eigengenes

# Compute module-trait correlations
module_trait_correlation <- cor(MEs, as.numeric(traits$KO_WT_status), use = "p")
module_trait_correlation_matrix <- as.matrix(module_trait_correlation)
colnames(module_trait_correlation_matrix) <- "KO_vs_WT"
rownames(module_trait_correlation_matrix) <- gsub("^ME", "", rownames(module_trait_correlation_matrix))

module_trait_pvalues <- corPvalueStudent(module_trait_correlation_matrix, nrow(counts_matrix_normalized))
order_idx <- order(abs(module_trait_correlation_matrix[, "KO_vs_WT"]), decreasing = TRUE)
sorted_correlation_matrix <- module_trait_correlation_matrix[order_idx, , drop = FALSE]
sorted_pvalues <- module_trait_pvalues[order_idx, , drop = FALSE]
text_matrix <- paste0(signif(sorted_correlation_matrix, 2), ifelse(sorted_pvalues < 0.05, "*", ""))
dim(text_matrix) <- dim(sorted_correlation_matrix)



labeledHeatmap(Matrix = sorted_correlation_matrix,
               xLabels = "KO_vs_WT",
               yLabels = rownames(sorted_correlation_matrix),
               ySymbols = rownames(sorted_correlation_matrix),
               colorLabels = FALSE,
               colors = blueWhiteRed(100),
               textMatrix = text_matrix,
               setStdMargins = TRUE,   # cleaner layout
               cex.text = 1,           # slightly larger, clearer labels
               cex.lab = 1.2,          # axis label size
               main = "",              # remove title for a clean top
               xLabelsAngle = 90,      # slight tilt for x axis
               zlim = c(-1, 1))        # standardize color scaling



# --- Create sample information for additional modeling ---

samples <- rownames(counts_matrix_normalized)
sample_info <- data.frame(
  Sample = samples,
  Status = ifelse(grepl("KO", samples), "KO", "WT"),
  CellType = case_when(
    grepl("^o", samples) ~ "Osteogenic",
    grepl("^ch", samples) ~ "Chondrogenic",
    grepl("^(WT|KO)", samples) ~ "Undifferentiated",
    TRUE ~ "Unknown"
  )
)
sample_info$Status <- factor(sample_info$Status, levels = c("WT", "KO"))
sample_info$CellType <- factor(sample_info$CellType, levels = c("Undifferentiated", "Osteogenic", "Chondrogenic"))
model_matrix <- model.matrix(~ Status + CellType, sample_info)[, -1]

# Module-trait correlation with the full model matrix
module_trait_correlation <- cor(MEs, model_matrix, use = "p")
module_trait_pvalues <- corPvalueStudent(module_trait_correlation, nSamples = nrow(counts_matrix_normalized))
text_matrix <- paste0(signif(module_trait_correlation, 2), "/n(",
                      signif(module_trait_pvalues, 1), ")")

col_fun <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
ht <- Heatmap(module_trait_correlation,
              name = "Correlation", col = col_fun,
              cluster_rows = TRUE, cluster_columns = TRUE,
              column_title = "Traits", row_title = "Modules",
              show_row_names = TRUE, show_column_names = TRUE,
              cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(sprintf("%.2f", module_trait_correlation[i, j]), x, y, gp = gpar(fontsize = 8))
              })
draw(ht, heatmap_legend_side = "right")

# (Optional) Plot module eigengene expression
pheatmap(as.matrix(MEs),
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Module Eigengene Expression")

# --- Calculate Adjacency & TOM using module eigengenes ---

MEs_clean <- MEs
if (any(sapply(MEs_clean, is.factor))) {
  MEs_clean <- as.data.frame(lapply(MEs_clean, function(x) as.numeric(as.character(x))))
}
beta <- 6
adjacency_matrix <- abs(cor(MEs_clean, use = "p"))^beta
TOM <- TOMsimilarity(adjacency_matrix)
colnames(TOM) <- colnames(MEs_clean)
rownames(TOM) <- colnames(MEs_clean)

# Prepare annotation for TOM heatmap
module_colors_modules <- gsub("^ME", "", colnames(MEs_clean))
annotation_row <- data.frame(Color = factor(module_colors_modules))
rownames(annotation_row) <- rownames(TOM)
annotation_col <- data.frame(Color = factor(module_colors_modules))
rownames(annotation_col) <- colnames(TOM)
annotation_colors <- list(Color = setNames(module_colors_modules, module_colors_modules))
custom_gradient <- colorRampPalette(c("#2E86C1", "white", "#CB4335"))(10)

pheatmap(TOM,
         color = custom_gradient,
         annotation_row = annotation_row,
         annotation_col = annotation_col,
         annotation_colors = annotation_colors,
         show_rownames = TRUE,
         show_colnames = TRUE,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         legend = TRUE,
         main = "Eigengene adjaceny",
         annotation_legend = FALSE,
         fontsize = 10)

# --- Map gene IDs for conversion (Ensembl -> HGNC) using biomaRt ---

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
conversion_table <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = gene_ids,  # original gene IDs from file
  mart = ensembl
)
conversion_table <- conversion_table[conversion_table$hgnc_symbol != "", ]
ensembl_to_symbol <- setNames(conversion_table$hgnc_symbol, conversion_table$ensembl_gene_id)
symbol_to_ensembl <- setNames(conversion_table$ensembl_gene_id, conversion_table$hgnc_symbol)

# --- Target Gene Analysis (Module Membership) ---

target_symbols <- c("FLI1", "WNT7A")

# Directly use the hard assignment from net$colors (recommended approach)
# Convert target gene symbols to Ensembl IDs first.
target_ensembl_ids <- symbol_to_ensembl[target_symbols]

# Use the hard assignment from blockwiseModules.
target_genes_net_assignment <- module_colors[target_ensembl_ids]
# Sanity check: print a data frame of target genes with their modules.
print(data.frame(Gene = target_symbols, Module = target_genes_net_assignment))

# (Optional) If you still want to compute module membership (MM) based on correlation for further insight,
# you can compute it here. However, for reporting the module assignment, we use the hard assignment.
module_membership <- cor(counts_matrix_normalized, MEs, use = "p")
target_genes_mm <- module_membership[target_ensembl_ids, ]
print("Module membership (MM) scores:")
print(target_genes_mm)

# --- Plot Target Gene Module Membership using Hard Assignment ---

# Prepare a data frame with gene names, MM (using the maximum correlation per gene) and the hard module assignment.
# We use the maximum MM from the computed correlation (which may differ slightly),
# but we annotate the bars using the hard assignment.
target_genes_mm_max <- apply(target_genes_mm, 1, max)
df_plot <- data.frame(
  Gene = target_symbols,
  Target_Gene_MM = target_genes_mm_max,
  Module = target_genes_net_assignment  # use the hard module assignment here
)

print("Data for target gene barplot:")
print(df_plot)

# Plot the target gene module membership with annotation of the module (color)
# For consistency, we use the module color from the hard assignment.
ggplot(df_plot, aes(x = Gene, y = Target_Gene_MM, fill = Module)) +
  geom_bar(stat = "identity", width = 0.7, color = "black", linewidth = 1) +  # <-- bar border here
  geom_text(aes(label = round(Target_Gene_MM, 2)),
            vjust = -0.5, size = 5, color = "black", fontface = "bold") +
  labs(title = "MM (FLI1, WNT7A)",
       y = "Maximum MM") +  # removed x title
  theme_classic() +
  theme(
    axis.text.x = element_text(hjust = 0.5, size = 18, face = "bold.italic"),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 16, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    legend.position = "none",
    axis.line = element_line(color = "black", linewidth = 1),
    panel.grid = element_blank(),
    panel.border = element_blank()
  ) +
  scale_fill_manual(values = setNames(unique(df_plot$Module), unique(df_plot$Module))) +
  ylim(0, 1)


# --- GO Enrichment Analysis for the "lightgreen" Module ---

# Extract the gene IDs (Ensembl) that are assigned to the lightgreen module.
genes_in_module <- names(module_colors)[module_colors == "lightgreen"]
# Check how many genes are in the "lightgreen" module
cat("Number of genes in the lightgreen module:", length(genes_in_module), "/n")
head(genes_in_module)

# Load the packages (if not already loaded)
library(clusterProfiler)
library(org.Hs.eg.db)

# Query biomaRt to convert Ensembl IDs to Entrez Gene IDs
gene_ids_conv <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id"),
                       filters = "ensembl_gene_id",
                       values = genes_in_module,
                       mart = ensembl)

# Clean up the Entrez Gene IDs (remove any NA values)
entrez_gene_ids <- na.omit(gene_ids_conv$entrezgene_id)
cat("First few Entrez Gene IDs:/n")
head(entrez_gene_ids)

# Perform GO term enrichment analysis with Entrez Gene IDs
go_results_BP_lightgreen <- enrichGO(gene = entrez_gene_ids,
                                     OrgDb = org.Hs.eg.db,
                                     ont = "BP",       # Change to "MF" or "CC" as needed
                                     pAdjustMethod = "BH",
                                     qvalueCutoff = 0.05,
                                     readable = TRUE)


library(dplyr)
library(ggplot2)
library(clusterProfiler)

df <- go_results_BP_lightgreen@result %>%  # extrahiere den DataFrame aus enrichResult-Objekt
  as_tibble() %>%
  arrange(p.adjust) %>%
  dplyr::slice(1:10) %>%
  mutate(
    Description = factor(Description, levels = rev(Description)),
    neg_log10_p = -log10(p.adjust)
  )


ggplot(df, aes(x = GeneRatio, y = Description, fill = neg_log10_p, size = Count)) +
  geom_point(shape = 21, color = "black", stroke = 0.5) +  # shape 21 = circle with border
  scale_fill_gradient(low = "lightgreen", high = "darkgreen", name = "-log10(p.adjust)") +
  scale_size(range = c(3, 8)) +
  labs(
    title = "GO-BP (lightgreen)",
    x = "Gene Ratio", y = NULL
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    axis.title.x = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14, face = "italic"),
    axis.line = element_line(color = "black", linewidth = 1),
    axis.ticks = element_line(color = "black", linewidth = 1),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    panel.grid = element_blank()
  )

# Extract result table
go_df <- as.data.frame(go_results_BP_lightgreen)

# Look for the term of interest
target_term <- "embryonic limb morphogenesis"

# Subset to the row for this GO term
embryonic_limb_row <- go_df[go_df$Description == target_term, ]

# Extract the gene symbols (semicolon separated string → split into vector)
if (nrow(embryonic_limb_row) > 0) {
  genes_embryonic_limb <- unlist(strsplit(embryonic_limb_row$geneID, "/"))
  cat("Genes enriched in", target_term, ":\n")
  print(genes_embryonic_limb)
} else {
  cat("GO term not found in enrichment results.\n")
}


library(ggplot2)
library(reshape2)
library(circlize)

# module_membership: matrix (genes × MEs)
kME_lightgreen <- module_membership[, "MElightgreen"]

n_top <- 15

top_hubs <- names(sort(abs(kME_lightgreen), decreasing = TRUE))[1:n_top]

mean_expr <- colMeans(counts_matrix_normalized[, top_hubs])
# re‐order or further filter if you like:
top_hubs <- names(sort(mean_expr, decreasing = TRUE))[1:n_top]


expr_sub <- counts_matrix_normalized[, top_hubs]
corr_sub  <- cor(expr_sub, use = "pairwise.complete.obs")

gene_syms_sub <- sapply(top_hubs, function(id) {
  sym <- ensembl_to_symbol[id]
  if (is.na(sym) || sym=="") id else sym
}, USE.NAMES = FALSE)
rownames(corr_sub) <- colnames(corr_sub) <- gene_syms_sub

library(pheatmap)
library(circlize)

# 1) Your precomputed hub‐gene correlation matrix, with HGNC symbols as row/colnames
#    → corr_sub

# 2) Define a Nature‐style palette (navy → white → goldenrod3)
my_cols <- colorRampPalette(c("navy", "white", "firebrick3"))(100)

# 3) (Optional) If you want a simple annotation of “high” vs “low” kME on the rows:
#    row_annotation <- data.frame(
#      kME = ifelse(abs(kME_lightgreen[top_hubs]) > 0.8, "hub", "non-hub")
#    )
#    rownames(row_annotation) <- rownames(corr_sub)
#    ann_colors <- list(kME = c(hub = "firebrick3", `non-hub` = "grey80"))

# --- 0. Manual mapping for missing / lncRNA genes ---
manual_symbols <- c(
  "ENSG00000240801" = "Lnc-INS-IGF2",        # Lnc-INS-IGF2-1 alias 
  "ENSG00000284779" = "ENSG00000284779", # uncharacterized / keep Ensembl as label (or change to a nicer alias)
  "ENSG00000291065" = "LOC441666",      # LOC441666
  "ENSG00000258973" = "Lnc-ADCY4-1"     # Lnc-ADCY4-1
)
# ensembl_to_symbol: a vector mapping Ensembl IDs -> HGNC symbols from previous annotation
# If not defined yet, you can create it with mapIds(...) like in your rescue plot

# --- 1. Replace Ensembl IDs with symbols / manual aliases ---
gene_syms_sub <- sapply(top_hubs, function(id) {
  sym <- ensembl_to_symbol[id]             # lookup from your annotation
  if (is.na(sym) || sym == "") {
    # fallback to manual mapping
    if (id %in% names(manual_symbols)) manual_symbols[id] else id
  } else {
    sym
  }
}, USE.NAMES = FALSE)

# --- 2. Update row and column names of your correlation matrix ---
rownames(corr_sub) <- colnames(corr_sub) <- gene_syms_sub

# --- 3. Plot as before ---
library(pheatmap)
my_cols <- colorRampPalette(c("navy", "white", "firebrick3"))(100)

pheatmap(
  mat               = corr_sub,
  color             = my_cols,
  name              = "Pearson r",
  cluster_rows      = TRUE,
  cluster_cols      = TRUE,
  show_rownames     = TRUE,
  show_colnames     = TRUE,
  fontsize          = 12,
  fontsize_row      = 12,
  fontsize_col      = 12,
  cellwidth         = 12,
  cellheight        = 12,
  border_color      = "white",
  legend_breaks     = c(-1, 0, 1),
  legend_labels     = c("-1", "0", "+1"),
  main              = "Hub genes (top 15)"
)

library(ComplexHeatmap)
ht <- Heatmap(
  corr_sub,
  name = "Pearson r",
  col  = colorRamp2(c(-1, 0, 1), c("navy", "white", "firebrick3")),
  row_names_gp = gpar(fontsize = 14, fontface = "bold.italic"),
  column_names_gp = gpar(fontsize = 14, fontface = "bold.italic"),
  rect_gp = gpar(col = "black", lwd = 0.5),
  column_title = "Top 15 hub genes",
  column_title_gp = gpar(fontsize = 18, fontface = "bold"),
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 14, fontface = "bold"),
    labels_gp = gpar(fontsize = 12)
  )
)



# Draw the heatmap first
ht_list <- draw(ht)

# Extract the row and column order
rows_in_heatmap <- row_order(ht_list)
cols_in_heatmap <- column_order(ht_list)

print(rows_in_heatmap)
print(cols_in_heatmap)

rownames(corr_sub)[rows_in_heatmap]
colnames(corr_sub)[cols_in_heatmap]


# Get expression vectors
fli1_expr <- counts_matrix_normalized[, names(ensembl_to_symbol[ensembl_to_symbol == "FLI1"])]
wnt7a_expr <- counts_matrix_normalized[, names(ensembl_to_symbol[ensembl_to_symbol == "WNT7A"])]

# Calculate Pearson correlation
cor_test <- cor.test(fli1_expr, wnt7a_expr, method = "pearson", use = "pairwise.complete.obs")

cor_test

df <- data.frame(
  FLI1 = fli1_expr,
  WNT7A = wnt7a_expr
)

ggplot(df, aes(x = FLI1, y = WNT7A)) +
  geom_point(color = "firebrick3", alpha = 0.7, size = 3, stroke = 0.5) +
  geom_smooth(method = "lm", se = FALSE, color = "navy", size = 1) +
  labs(
    title = "FLI1 vs WNT7A",
    subtitle = paste0("r = ",
                      round(cor_test$estimate, 2),
                      ", p = ", signif(cor_test$p.value, 3)),
    x = "FLI1",
    y = "WNT7A"
  ) +
  theme_void(base_family = "sans") +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5, margin = margin(b = 10)),
    plot.subtitle = element_text(size = 14, face = "plain", hjust = 0.5, margin = margin(b = 10)), # plain subtitle
    plot.margin = margin(30, 20, 20, 20),
    axis.title.y = element_text(
      angle = 90,              
      size = 16,               
      face = "bold.italic",
      vjust = 5
    ),
    axis.title.x = element_text(size = 16, face = "bold.italic", vjust = - 3),
    axis.line = element_line(color = "black", linewidth = 1),
    axis.ticks = element_line(color = "black", linewidth = 1)
  )



library(WGCNA)
library(igraph)
library(ggraph)
library(ggplot2)
library(ggrepel)

# --- Extract lightgreen genes ---
genes_lightgreen <- names(module_colors)[module_colors == "lightgreen"]

# Expression restricted to this module
expr_lightgreen <- counts_matrix_normalized[, genes_lightgreen]

# --- Adjacency and TOM restricted to lightgreen genes ---
adjacency_lightgreen <- adjacency(expr_lightgreen, power = chosen_power, type = "unsigned")
TOM_lightgreen <- TOMsimilarity(adjacency_lightgreen)

colnames(TOM_lightgreen) <- rownames(TOM_lightgreen) <- genes_lightgreen

# --- Build igraph object from adjacency/TOM ---
# You can use adjacency or TOM (TOM is usually better for network structure)
threshold <- 0.25 # keep only stronger edges
adj_thr <- TOM_lightgreen
adj_thr[adj_thr < threshold] <- 0

g_lightgreen <- graph.adjacency(adj_thr, mode = "undirected", weighted = TRUE, diag = FALSE)

# Remove isolated nodes
g_lightgreen <- delete.vertices(g_lightgreen, which(degree(g_lightgreen) == 0))
cat("Vertices remaining in graph:", vcount(g_lightgreen), "\n")

# --- Add vertex attributes ---
V(g_lightgreen)$ensembl <- as.character(V(g_lightgreen)$name)
V(g_lightgreen)$symbol <- ensembl_to_symbol[V(g_lightgreen)$ensembl]
V(g_lightgreen)$symbol[is.na(V(g_lightgreen)$symbol)] <- V(g_lightgreen)$ensembl[is.na(V(g_lightgreen)$symbol)]

# --- Highlight FLI1 (and optionally WNT7A) ---
highlight_genes <- c("FLI1", "WNT7A")
V(g_lightgreen)$highlight <- V(g_lightgreen)$symbol %in% highlight_genes

# --- Define GO term genes ---
go_genes <- c("PITX1", "GRHL2", "TBX2", "TFAP2A", "WNT7A")
go_genes <- intersect(go_genes, V(g_lightgreen)$symbol)  # keep only genes in graph

library(WGCNA)
library(igraph)
library(ggraph)
library(ggplot2)
library(ggrepel)

# --- Extract lightgreen genes ---
genes_lightgreen <- names(module_colors)[module_colors == "lightgreen"]

# Expression restricted to this module
expr_lightgreen <- counts_matrix_normalized[, genes_lightgreen]

# --- Adjacency and TOM restricted to lightgreen genes ---
adjacency_lightgreen <- adjacency(expr_lightgreen, power = chosen_power, type = "unsigned")
TOM_lightgreen <- TOMsimilarity(adjacency_lightgreen)

colnames(TOM_lightgreen) <- rownames(TOM_lightgreen) <- genes_lightgreen

# --- Build igraph object from adjacency/TOM ---
# You can use adjacency or TOM (TOM is usually better for network structure)
threshold <- 0.25  # keep only stronger edges
adj_thr <- TOM_lightgreen
adj_thr[adj_thr < threshold] <- 0

g_lightgreen <- graph.adjacency(adj_thr, mode = "undirected", weighted = TRUE, diag = FALSE)

# Remove isolated nodes
g_lightgreen <- delete.vertices(g_lightgreen, which(degree(g_lightgreen) == 0))
cat("Vertices remaining in graph:", vcount(g_lightgreen), "\n")

# --- Add vertex attributes ---
V(g_lightgreen)$ensembl <- as.character(V(g_lightgreen)$name)
V(g_lightgreen)$symbol <- ensembl_to_symbol[V(g_lightgreen)$ensembl]
V(g_lightgreen)$symbol[is.na(V(g_lightgreen)$symbol)] <- V(g_lightgreen)$ensembl[is.na(V(g_lightgreen)$symbol)]

# --- Highlight FLI1 (and optionally WNT7A) ---
highlight_genes <- c("FLI1", "WNT7A")
V(g_lightgreen)$highlight <- V(g_lightgreen)$symbol %in% highlight_genes

V(g_lightgreen)$color <- ifelse(V(g_lightgreen)$highlight, "red", "lightgreen")
V(g_lightgreen)$frame.color <- ifelse(V(g_lightgreen)$highlight, "black", "gray70")
V(g_lightgreen)$frame.width <- ifelse(V(g_lightgreen)$highlight, 2.5, 1)
V(g_lightgreen)$label <- V(g_lightgreen)$symbol

# --- Nature-style ggraph plot ---
ggraph(g_lightgreen, layout = "sphere") +
  geom_edge_link(aes(width = weight), alpha = 0.3, color = "gray60") +
  geom_node_point(aes(color = I(color)),
                  size = 5,
                  stroke = 0.6,
                  shape = 21,
                  color = "black",
                  fill = V(g_lightgreen)$color) +
  geom_node_text(
    aes(label = label),
    repel = TRUE,
    size = 5,                 
    fontface = "bold.italic", 
    box.padding = unit(0.4, "lines"),
    point.padding = unit(0.3, "lines"),
    segment.size = 0.2,
    family = "sans",
    max.overlaps = 50,
    force = 1.2
  ) +
  scale_edge_width(range = c(0.2, 1.5)) +
  ggtitle("Lightgreen network (TOM > 0.25)") +
  theme_void(base_family = "sans") +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5, margin = margin(b = 10)),
    plot.margin = margin(30, 20, 20, 20),
    legend.position = "none"
  )



# Script A: Export lightgreen module genes to Excel

library(writexl)
library(biomaRt)       # for mapping Ensembl -> HGNC (optional)
library(dplyr)



# Extract lightgreen genes:
genes_lightgreen <- names(module_colors)[module_colors == "lightgreen"]
cat("Lightgreen genes found:", length(genes_lightgreen), "\n")

# Map to HGNC symbols via biomaRt (optional)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
conversion_table <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = genes_lightgreen,
  mart = ensembl
)

# Keep unique mapping (some Ensembl->multiple symbols can exist; we take first)
conversion_table <- conversion_table %>%
  dplyr::group_by(ensembl_gene_id) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup()

library(dplyr)
library(writexl)

# --- Compute TOM-based connectivity ---
# For each gene, sum all TOM connections above a threshold (excluding self-connections)
tom_threshold <- 0.25
diag(TOM_lightgreen) <- 0  # exclude self-connections
gene_tom_score <- rowSums(TOM_lightgreen > tom_threshold)  # number of strong connections

# --- Define hub genes ---
hub_genes <- names(gene_tom_score[gene_tom_score > 0])  # genes with at least one TOM connection > threshold

# --- Build final table ---
out_tbl <- data.frame(
  ensembl_gene_id = genes_lightgreen,
  stringsAsFactors = FALSE
) %>%
  left_join(conversion_table, by = c("ensembl_gene_id" = "ensembl_gene_id")) %>%
  mutate(
    hgnc_symbol = ifelse(is.na(hgnc_symbol) | hgnc_symbol == "",
                         ensembl_gene_id,
                         hgnc_symbol),
    is_hub = ensembl_gene_id %in% hub_genes
  )

# --- Write to Excel ---
write_xlsx(list("lightgreen_genes" = out_tbl),
           path = "lightgreen_genes_with_hub_TOM.xlsx")
cat("Written lightgreen_genes_with_hub_TOM.xlsx (", nrow(out_tbl), "rows )\n", sep = "")
