# Load libraries
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(org.Hs.eg.db)
library(AnnotationDbi)

# --- Step 1: Load data ---
file_path <- "...Counts_ThielRNA.xlsx"
raw_counts <- read_excel(file_path, sheet = 1)

# --- Step 2: Prepare counts (select only desired samples) ---
selected_samples <- c("WT_1", "WT_2", "WT_3", "KO_9", "KO_75", "KO_46")
counts_matrix <- as.matrix(raw_counts[, selected_samples])
rownames(counts_matrix) <- raw_counts$ID

# --- Step 3: Map Ensembl IDs to gene symbols ---
gene_symbols <- mapIds(org.Hs.eg.db,
                       keys = rownames(counts_matrix),
                       column = "SYMBOL",
                       keytype = "ENSEMBL",
                       multiVals = "first")

# Install if needed
if (!requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE)) {
  BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
}

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# Sum exon lengths per gene
exons_by_gene <- exonsBy(txdb, by="gene")
# Sum exon widths per gene → gives a named vector of lengths
gene_lengths_bp <- sapply(exons_by_gene, function(x) sum(width(reduce(x))))
gene_lengths_kb <- gene_lengths_bp / 1000

# Map Entrez IDs to your Ensembl IDs
library(org.Hs.eg.db)
# Remove version from Ensembl IDs
ensembl_ids <- sub("\\..*", "", rownames(counts_matrix))

# Map Ensembl IDs to Entrez IDs
entrez_ids <- mapIds(org.Hs.eg.db, keys=ensembl_ids,
                     column="ENTREZID", keytype="ENSEMBL", multiVals="first")

# Align lengths to counts
gene_lengths_kb <- gene_lengths_kb[entrez_ids]
valid_genes <- !is.na(gene_lengths_kb)
counts_matrix <- counts_matrix[valid_genes, ]
gene_lengths_kb <- gene_lengths_kb[valid_genes]

# TPM function
tpm <- function(counts, lengths_kb) {
  rate <- counts / lengths_kb
  t(t(rate) / colSums(rate) * 1e6)
}

tpm_matrix <- tpm(counts_matrix, gene_lengths_kb)


# --- Step 6: Prepare data for plotting ---
plot_data <- as.data.frame(tpm_matrix) %>%
  mutate(Ensembl_ID = rownames(tpm_matrix),
         Gene_Symbol = gene_symbols[rownames(tpm_matrix)]) %>%
  pivot_longer(cols = selected_samples, names_to = "Sample", values_to = "TPM") %>%
  mutate(Condition = ifelse(grepl("WT", Sample), "WT", "KO"))

# --- Step 7: Summarize for mean ± SD per condition ---
genes_of_interest <- c("PITX1", "GRHL2", "TBX2", "TFAP2A", "WNT7A")  

plot_summary <- plot_data %>%
  filter(Gene_Symbol %in% genes_of_interest) %>%
  group_by(Gene_Symbol, Condition) %>%
  summarise(mean_TPM = mean(TPM),
            sd_TPM = sd(TPM),
            .groups = "drop")


# --- Step 8: Plot merged bars ---
# --- Styled plot similar to module membership plot ---
# --- Styled plot with requested changes ---
# --- p1: error bars only above ---
p1 <- ggplot(plot_summary, aes(x = Condition, y = mean_TPM, fill = Condition)) +
  geom_bar(stat = "identity", width = 0.7, color = "black", linewidth = 1) +
  geom_errorbar(aes(ymin = mean_TPM, ymax = mean_TPM + sd_TPM),
                width = 0.2, linewidth = 1) +
  # >>> ADD THIS <<<
  geom_jitter(
    data = plot_data %>% filter(Gene_Symbol %in% genes_of_interest),
    aes(x = Condition, y = TPM, fill = Condition),
    width = 0.15,
    size = 2.5,
    shape = 21,
    color = "black",
    inherit.aes = FALSE
  ) +
  facet_wrap(~Gene_Symbol, scales = "free_y", strip.position = "top") +
  labs(y = "TPM", x = "") +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold.italic", size = 14),
    axis.text.x = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14, face = "bold"),
    legend.position = "none",
    axis.line = element_line(color = "black", linewidth = 1),
    panel.grid = element_blank()
  ) +
  scale_fill_manual(values = c("WT" = "lightgreen", "KO" = "#f66666"))


# Load libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(org.Hs.eg.db)
library(AnnotationDbi)

# --- Step 1: Load .tab data ---
file_path <- "...ThielRNA.tab"

raw_counts <- read.delim(file_path, header = TRUE, stringsAsFactors = FALSE)

# --- Step 2: Select only desired samples ---
# KO samples: ThielRNA_KO1/3/4
# LiCl samples: ThielRNA_LiCl1/2/3
selected_samples <- c("ThielRNA_KO1.bam", "ThielRNA_KO3.bam", "ThielRNA_KO4.bam",
                      "ThielRNA_LiCl1.bam", "ThielRNA_LiCl2.bam", "ThielRNA_LiCl3.bam")

counts_matrix <- as.matrix(raw_counts[, selected_samples])
rownames(counts_matrix) <- raw_counts$Geneid

# --- Step 3: Map Ensembl IDs to gene symbols ---
gene_symbols <- mapIds(org.Hs.eg.db,
                       keys = rownames(counts_matrix),
                       column = "SYMBOL",
                       keytype = "ENSEMBL",
                       multiVals = "first")

# Install if needed
if (!requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE)) {
  BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
}

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# Sum exon lengths per gene
exons_by_gene <- exonsBy(txdb, by="gene")
# Sum exon widths per gene → gives a named vector of lengths
gene_lengths_bp <- sapply(exons_by_gene, function(x) sum(width(reduce(x))))
gene_lengths_kb <- gene_lengths_bp / 1000


# Map Entrez IDs to your Ensembl IDs
library(org.Hs.eg.db)
# Remove version from Ensembl IDs
ensembl_ids <- sub("\\..*", "", rownames(counts_matrix))

# Map Ensembl IDs to Entrez IDs
entrez_ids <- mapIds(org.Hs.eg.db, keys=ensembl_ids,
                     column="ENTREZID", keytype="ENSEMBL", multiVals="first")

# Align lengths to counts
gene_lengths_kb <- gene_lengths_kb[entrez_ids]
valid_genes <- !is.na(gene_lengths_kb)
counts_matrix <- counts_matrix[valid_genes, ]
gene_lengths_kb <- gene_lengths_kb[valid_genes]

# TPM function
tpm <- function(counts, lengths_kb) {
  rate <- counts / lengths_kb
  t(t(rate) / colSums(rate) * 1e6)
}

tpm_matrix <- tpm(counts_matrix, gene_lengths_kb)


# --- Step 6: Prepare data for plotting ---
plot_data <- as.data.frame(tpm_matrix) %>%
  mutate(Ensembl_ID = rownames(tpm_matrix),
         Gene_Symbol = gene_symbols[rownames(tpm_matrix)]) %>%
  pivot_longer(cols = selected_samples, names_to = "Sample", values_to = "TPM") %>%
  mutate(Condition = ifelse(grepl("KO", Sample), "KO", "LiCl"))



plot_summary <- plot_data %>%
  filter(Gene_Symbol %in% genes_of_interest) %>%
  group_by(Gene_Symbol, Condition) %>%
  summarise(mean_TPM = mean(TPM),
            sd_TPM = sd(TPM),
            .groups = "drop")


# --- p2: error bars only above ---
p2 <- ggplot(plot_summary, aes(x = Condition, y = mean_TPM, fill = Condition)) +
  geom_bar(stat = "identity", width = 0.7, color = "black", linewidth = 1) +
  geom_errorbar(aes(ymin = mean_TPM, ymax = mean_TPM + sd_TPM),
                width = 0.2, linewidth = 1) +
  # >>> ADD THIS <<<
  geom_jitter(
    data = plot_data %>% filter(Gene_Symbol %in% genes_of_interest),
    aes(x = Condition, y = TPM, fill = Condition),
    width = 0.15,
    size = 2.5,
    shape = 21,
    color = "black",
    inherit.aes = FALSE
  ) +
  facet_wrap(~Gene_Symbol, scales = "free_y", strip.position = "top") +
  labs(y = "TPM", x = "") +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold.italic", size = 14),
    axis.text.x = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14, face = "bold"),
    legend.position = "none",
    axis.line = element_line(color = "black", linewidth = 1),
    panel.grid = element_blank()
  ) +
  scale_fill_manual(values = c("KO" = "#f66666", "LiCl" = "lightgreen"))


library(patchwork)

# --- add per-plot titles and styling ---
p1_titled <- p1 +
  labs(title = "KO vs. WT") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18)
  )

p2_titled <- p2 +
  labs(title = "KO2 vs. KO2 (LiCl)") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18)
  )

# --- keep facets stacked so each facet shows its own x-axis ---
p1_ncol1 <- p1_titled + facet_wrap(~Gene_Symbol, scales = "free_y", ncol = 1) +
  geom_hline(yintercept = 0, color = "black", linewidth = 1)
p2_ncol1 <- p2_titled + facet_wrap(~Gene_Symbol, scales = "free_y", ncol = 1) +
  geom_hline(yintercept = 0, color = "black", linewidth = 1)

# --- combine and add the big centered figure title ---
combined <- p1_ncol1 - p2_ncol1 +
  plot_annotation(
    title = "Embryonic limb morphogenesis rescue",
    theme = theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 20)
    )
  )

# print
combined


library(patchwork)

# --- adjust facets to be in one row per plot ---
p1_nrow1 <- p1_titled + facet_wrap(~Gene_Symbol, scales = "free_y", nrow = 1)
p2_nrow1 <- p2_titled + facet_wrap(~Gene_Symbol, scales = "free_y", nrow = 1) 
# --- combine vertically ---
combined <- p1_nrow1 / p2_nrow1 +
  plot_annotation(
    title = "Embryonic limb morphogenesis gene expression rescue",
    theme = theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 20)
    )
  )

# print
combined

# Load your DESeq2 results
res_KO_WT <- read_excel("...differentialExpression_Thiel_naiv.xlsx")
res_KO_LiCl <- read_excel("...DESeq2_KO_vs_LiCl_results.xlsx")



# ---------- significance annotation (add after you read res_KO_WT / res_KO_LiCl) ----------
library(dplyr)

# helper: map padj -> stars
padj_to_stars <- function(p) {
  dplyr::case_when(
    is.na(p)            ~ NA_character_,
    p < 0.001           ~ "***",
    p < 0.01            ~ "**",
    p < 0.05            ~ "*",
    TRUE                ~ NA_character_
  )
}

# Prepare annotation for p1 (KO vs WT)
ann_p1 <- plot_summary %>%
  # keep only genes-of-interest (plot_summary already filtered earlier, but safe)
  filter(Gene_Symbol %in% genes_of_interest) %>%
  # take the KO rows (we will place the star above the KO bar)
  filter(Condition == "KO") %>%
  left_join(
    res_KO_WT %>% dplyr::select(gene_symbol, padj),
    by = c("Gene_Symbol" = "gene_symbol")
  ) %>%
  mutate(
    sig = padj_to_stars(padj),
    # place the star above the errorbar: mean + sd + small offset
    y = mean_TPM + sd_TPM + pmax(sd_TPM * 0.4, 0.05)
  ) %>%
  filter(!is.na(sig))

# Prepare annotation for p2 (KO vs LiCl)
ann_p2 <- plot_summary %>%
  filter(Gene_Symbol %in% genes_of_interest) %>%
  filter(Condition == "KO") %>%
  left_join(
    res_KO_LiCl %>% dplyr::select(gene_symbol, padj),
    by = c("Gene_Symbol" = "gene_symbol")
  ) %>%
  mutate(
    sig = padj_to_stars(padj),
    y = mean_TPM + sd_TPM + pmax(sd_TPM * 0.4, 0.05)
  ) %>%
  filter(!is.na(sig))

# Add the annotation to existing ggplots (these objects exist in your script: p1 and p2)
p1 <- p1 +
  geom_text(
    data = ann_p1,
    aes(x = Condition, y = y, label = sig),
    inherit.aes = FALSE,
    size = 6,         # change size if you want larger/smaller stars
    vjust = 0         # place the text exactly at y (above bar)
  )

p2 <- p2 +
  geom_text(
    data = ann_p2,
    aes(x = Condition, y = y, label = sig),
    inherit.aes = FALSE,
    size = 6,
    vjust = 0
  )

library(patchwork)

# --- add per-plot titles and styling ---
p1_titled <- p1 +
  labs(title = "KO vs. WT") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20)
  )

p2_titled <- p2 +
  labs(title = "KO vs. KO (LiCl)") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20)
  )

# --- keep facets stacked so each facet shows its own x-axis ---
p1_ncol1 <- p1_titled + facet_wrap(~Gene_Symbol, scales = "free_y", ncol = 1) +
  geom_hline(yintercept = 0, color = "black", linewidth = 1)
p2_ncol1 <- p2_titled + facet_wrap(~Gene_Symbol, scales = "free_y", ncol = 1) +
  geom_hline(yintercept = 0, color = "black", linewidth = 1)

# --- combine and add the big centered figure title ---
combined <- p1_ncol1 - p2_ncol1 +
  plot_annotation(
    title = "Embryonic limb morphogenesis gene expression",
    theme = theme(
      plot.title = element_text(hjust = 0.5, face = "italic", size = 24)
    )
  )

# print
combined


library(patchwork)

# --- adjust facets to be in one row per plot ---
p1_nrow1 <- p1_titled + facet_wrap(~Gene_Symbol, scales = "free_y", nrow = 1)
p2_nrow1 <- p2_titled + facet_wrap(~Gene_Symbol, scales = "free_y", nrow = 1) 
# --- combine vertically ---
combined <- p1_nrow1 / p2_nrow1 +
  plot_annotation(
    title = "Embryonic limb morphogenesis gene expression",
    theme = theme(
      plot.title = element_text(hjust = 0.5, face = "italic", size = 24)
    )
  )

# print
combined

