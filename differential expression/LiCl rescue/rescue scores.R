# =============================
# LiCl Rescue Analysis (batch-corrected VST)
# =============================

# --- Libraries ---
library(DESeq2)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readxl)
library(limma)      # removeBatchEffect
library(scales)     # percent_format

# --- 1. Load counts ---
# Batch1: WT + KO2 (main experiment)
file_batch1 <- "...Counts_ThielRNA.xlsx"
batch1_df <- read_excel(file_batch1)
samples_batch1 <- c("WT_1","WT_2","WT_3","KO_75")  # WT + KO2
counts1 <- as.matrix(batch1_df[, samples_batch1])
rownames(counts1) <- batch1_df$ID

bio_map1 <- c("WT_1"="WT","WT_2"="WT","WT_3"="WT","KO_75"="KO2")
bio1 <- unname(bio_map1[colnames(counts1)])
treatment1 <- rep("none", ncol(counts1))
batch1 <- rep("Batch1", ncol(counts1))
coldata1 <- data.frame(row.names=colnames(counts1),
                       bio=bio1,
                       treatment=treatment1,
                       batch=batch1)

# Batch2: KO2 ± LiCl (rescue)
file_batch2 <- "...ThielRNA.tab"
batch2_df <- read.delim(file_batch2, header=TRUE, stringsAsFactors=FALSE)
counts2 <- as.matrix(batch2_df[, grep("ThielRNA", colnames(batch2_df))])
rownames(counts2) <- batch2_df$Geneid
batch2 <- rep("Batch2", ncol(counts2))
# Detect treatment from column names if possible
treatment2 <- ifelse(grepl("LiCl", colnames(counts2), ignore.case = TRUE), "LiCl", "none")
bio2 <- rep("KO2", ncol(counts2))
coldata2 <- data.frame(row.names=colnames(counts2),
                       bio=bio2,
                       treatment=treatment2,
                       batch=batch2)

# --- 2. Harmonize genes and combine ---
common_genes <- intersect(rownames(counts1), rownames(counts2))
counts1 <- counts1[common_genes, ]
counts2 <- counts2[common_genes, ]
combined_counts <- cbind(counts1, counts2)
combined_coldata <- rbind(coldata1, coldata2)
combined_coldata <- combined_coldata[colnames(combined_counts), ]
combined_coldata$batch <- factor(combined_coldata$batch)
combined_coldata$bio <- factor(combined_coldata$bio)
combined_coldata$treatment <- factor(combined_coldata$treatment)

# --- 3. DESeq2 for differential testing ---
dds <- DESeqDataSetFromMatrix(round(combined_counts),
                              colData=combined_coldata,
                              design=~ batch + bio + treatment)
dds <- DESeq(dds)

# --- 4. VST transformation ---
vsd <- vst(dds, blind=FALSE)

# --- 5. Batch correction (removeBatchEffect) ---
design_matrix <- model.matrix(~ bio + treatment, data=as.data.frame(colData(vsd)))
mat_vst_bc <- removeBatchEffect(assay(vsd), batch=combined_coldata$batch, design=design_matrix)

# --- 6. Gene symbols & genes of interest ---
gene_symbols <- mapIds(org.Hs.eg.db,
                       keys=rownames(combined_counts),
                       column="SYMBOL",
                       keytype="ENSEMBL",
                       multiVals="first")
genes_of_interest <- c("PITX1", "GRHL2", "TBX2", "TFAP2A", "WNT7A")
ensembl_selected <- names(gene_symbols)[gene_symbols %in% genes_of_interest]
gene_symbols_selected <- gene_symbols[ensembl_selected]

# --- 7. Long-format VST for plotting ---
gene_idx <- rownames(mat_vst_bc)[rownames(mat_vst_bc) %in% ensembl_selected]
vst_long <- as.data.frame(mat_vst_bc[gene_idx, ]) %>%
  tibble::rownames_to_column("Gene") %>%
  mutate(Gene = gene_symbols_selected[Gene]) %>%
  pivot_longer(-Gene, names_to="Sample", values_to="Expression") %>%
  left_join(as.data.frame(combined_coldata) %>% tibble::rownames_to_column("Sample"), by="Sample") %>%
  mutate(Group = case_when(
    bio == "WT" ~ "WT",
    bio == "KO2" & treatment == "LiCl" ~ "KO2+LiCl",
    bio == "KO2" & treatment != "LiCl" ~ "KO2"
  ),
  Group = factor(Group, levels=c("WT","KO2","KO2+LiCl"))
  )

# --- 8. Compute batch-corrected rescue index ---
rescue_df <- vst_long %>%
  group_by(Gene, Group) %>%
  summarise(mean_expr=mean(Expression), .groups="drop") %>%
  pivot_wider(names_from=Group, values_from=mean_expr) %>%
  mutate(rescue_index = (`KO2+LiCl` - KO2) / (WT - KO2))

# --- 9. Rescue plot ---
p_rescue <- ggplot(rescue_df, aes(x=Gene, y=rescue_index, fill=rescue_index)) +
  geom_col(fill="gray80", color="black", linewidth=1) +  # uniform gray fill
  geom_hline(yintercept=1, linetype="dashed", color="gray40", linewidth=1) +
  geom_hline(yintercept=0, linetype="dotted", color="gray60", linewidth=1) +
  coord_flip() +
  geom_text(aes(label = scales::percent(rescue_index, accuracy = 1)), 
            hjust = -0.1, size = 5, fontface="bold.italic") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  theme_classic(base_size=16) +
  theme(
    axis.text.y = element_text(face = "bold.italic", size=16),
    axis.text.x = element_text(size=14),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size=16),
    plot.title = element_text(size=18, face="bold"),
    plot.subtitle = element_text(size=16, face="italic"),
    legend.position = "none"
  ) +
  labs(y="Rescue index (%)", x="", title="LiCl treatment KO", subtitle = "Embryonic limb morphogenesis")

# --- 10. Print rescue plot ---
print(p_rescue)

# --- add manual mappings for lncRNAs / missing symbols ---
# adjust these names/aliases to taste
manual_symbols <- c(
  "ENSG00000240801" = "Lnc-INS-IGF2",        # Lnc-INS-IGF2-1 alias -> use "INS-IGF2" to match your genes_of_interest
  "ENSG00000284779" = "ENSG00000284779", # uncharacterized / keep Ensembl as label (or change to a nicer alias)
  "ENSG00000291065" = "LOC441666",      # LOC441666
  "ENSG00000258973" = "Lnc-ADCY4-1"     # Lnc-ADCY4-1
)

# gene_symbols was created earlier with mapIds(...)
# create a final symbol lookup that prefers org.Hs.eg.db, then manual, then Ensembl ID
gene_symbols_final <- gene_symbols                     # start with mapped symbols
# apply manual overrides / fills
for(id in names(manual_symbols)) {
  gene_symbols_final[id] <- manual_symbols[id]
}
# fallback: if still NA, fill with Ensembl ID itself
na_idx <- is.na(gene_symbols_final) | gene_symbols_final == ""
gene_symbols_final[na_idx] <- names(gene_symbols_final)[na_idx]

# --- genes of interest: keep your original desired set (symbols) ---
genes_of_interest <- c("IGF2", "INS-IGF2", "MARCHF3", "POPDC3",
                       "FLI1", "VEGFC", "PLD5", "SLC22A2",
                       "PCDHB2", "WNT7A", "PTGER3")

# --- choose Ensembl IDs that match either the symbol list OR are explicitly in manual_symbols ---
ensembl_selected <- names(gene_symbols_final)[
  gene_symbols_final %in% genes_of_interest | names(gene_symbols_final) %in% names(manual_symbols)
]
# (This includes Ensembl IDs whose assigned alias is in genes_of_interest,
#  and also ensures we include the manual ENSG IDs you added.)

# Subset the mapping to the selected set
gene_symbols_selected <- gene_symbols_final[ensembl_selected]

# --- 7. Long-format VST for plotting (use batch-corrected mat_vst_bc) ---
# note: mat_vst_bc is the removeBatchEffect-corrected VST matrix from earlier
gene_idx <- rownames(mat_vst_bc)[rownames(mat_vst_bc) %in% ensembl_selected]

vst_long <- as.data.frame(mat_vst_bc[gene_idx, , drop = FALSE]) %>%
  tibble::rownames_to_column("GeneID") %>%
  pivot_longer(cols = -GeneID, names_to = "Sample", values_to = "Expression") %>%
  # map GeneID -> human-friendly label (use gene_symbols_selected, fallback to GeneID)
  mutate(GeneSymbol = gene_symbols_selected[GeneID],
         Gene = ifelse(is.na(GeneSymbol) | GeneSymbol == "", GeneID, GeneSymbol)) %>%
  dplyr::select(-GeneSymbol) %>%
  left_join(as.data.frame(combined_coldata) %>% tibble::rownames_to_column("Sample"), by = "Sample") %>%
  mutate(Group = case_when(
    bio == "WT" ~ "WT",
    bio == "KO2" & treatment == "LiCl" ~ "KO2+LiCl",
    bio == "KO2" & treatment != "LiCl" ~ "KO2"
  ),
  Group = factor(Group, levels = c("WT", "KO2", "KO2+LiCl"))
  )

# --- 8. Compute batch-corrected rescue index ---
rescue_df <- vst_long %>%
  group_by(Gene, Group) %>%
  summarise(mean_expr = mean(Expression), .groups = "drop") %>%
  pivot_wider(names_from = Group, values_from = mean_expr) %>%
  # if any of the required columns are missing (e.g. no KO2+LiCl samples for a gene), keep NA safe
  mutate(rescue_index = (`KO2+LiCl` - KO2) / (WT - KO2))
# --- 1. Define the desired plotting order from your heatmap ---
heatmap_order <- c(
  "Lnc-INS-IGF2", "IGF2", "ENSG00000284779", "INS-IGF2", "MARCHF3",
  "POPDC3", "FLI1", "VEGFC", "PLD5", "SLC22A2",
  "PCDHB2", "WNT7A", "LOC441666", "PTGER3", "Lnc-ADCY4-1"
)

# --- 2. Ensure rescue_df$Gene is a factor in this order ---
rescue_df$Gene <- factor(rescue_df$Gene, levels = rev(heatmap_order))


# --- 3. Plot with uniform gray fill, respecting the heatmap order ---
p_rescue <- ggplot(rescue_df, aes(x = Gene, y = rescue_index)) +
  geom_col(fill = "gray80", color = "black", linewidth = 1) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40", linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "gray60", linewidth = 1) +
  coord_flip() +
  geom_text(aes(label = scales::percent(rescue_index, accuracy = 1)),
            hjust = -0.1, size = 5, fontface = "bold.italic") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  theme_classic(base_size = 16) +
  theme(
    axis.text.y = element_text(face = "bold.italic", size = 16),
    axis.text.x = element_text(size = 14),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 16),
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 16, face = "italic"),
    legend.position = "none"
  ) +
  labs(y = "Rescue index (%)", x = "", title = "LiCl treatment KO", subtitle = "Top 15 hub genes")

# --- 4. Print the plot ---
print(p_rescue)
