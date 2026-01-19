# === Libraries ===
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(readxl)
library(limma)    # for removeBatchEffect
# library(sva)    # optional: for ComBat if you prefer

# === Files (update paths if needed) ===
file1 <- "...ThielRNA.tab"
file2 <- "...Counts_ThielRNA.xlsx"

# === 1. Load first count table (rescue / LiCl experiment, which you've been calling "Batch2") ===
count_data1 <- read.table(file1, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
rownames(count_data1) <- count_data1$Geneid
counts1 <- count_data1[, -1, drop = FALSE]

# For counts1 the user told us: all samples are KO2; 3 untreated, 3 LiCl-treated.
# We'll assign sample names as they appear in the count table columns.
batch1 <- rep("Batch2", ncol(counts1))
# Create provisional colnames for these rescue samples if they are messy; we'll keep original names
condition1_raw <- colnames(counts1)   # these likely are e.g. ThielRNA_KO1.bam etc
# We'll set bio = "KO2" for all, and treatment depending on sample order (user said first 3 untreated, next 3 LiCl)
# If your file has explicit naming, adjust accordingly.
n1 <- ncol(counts1)
# blank-check: if not exactly 6, we still try to detect LiCl in names
if(n1 == 6) {
  treatment1 <- c(rep("none", 3), rep("LiCl", 3))
} else {
  # try to detect "LiCl" in column names; fallback: mark all 'none'
  treatment1 <- ifelse(grepl("LiCl", condition1_raw, ignore.case = TRUE), "LiCl", "none")
}
coldata1 <- data.frame(row.names = colnames(counts1),
                       bio = rep("KO2", ncol(counts1)),
                       treatment = treatment1,
                       batch = batch1,
                       stringsAsFactors = FALSE)

# === 2. Load second count table (main experiment, "Batch1") ===
count_data2 <- read_excel(file2)
count_data2 <- as.data.frame(count_data2)

# Keep only samples of interest (as in your original code)
samples_of_interest <- c("WT_1","WT_2","WT_3","KO_9","KO_75","KO_46")
counts2 <- count_data2[, c("ID", samples_of_interest)]
rownames(counts2) <- counts2$ID
counts2 <- counts2[, -1, drop = FALSE]

# For the main experiment:
# Map the KO sample identifiers to meaningful biological labels:
# KO_9  -> KO1
# KO_75 -> KO2  (this is the real bridge sample)
# KO_46 -> KO3
bio_map <- c("WT_1" = "WT", "WT_2" = "WT", "WT_3" = "WT",
             "KO_9" = "KO1", "KO_75" = "KO2", "KO_46" = "KO3")
bio2 <- unname(bio_map[colnames(counts2)])
batch2 <- rep("Batch1", ncol(counts2))
treatment2 <- rep("none", ncol(counts2))  # main experiment had no LiCl treatment
coldata2 <- data.frame(row.names = colnames(counts2),
                       bio = bio2,
                       treatment = treatment2,
                       batch = batch2,
                       stringsAsFactors = FALSE)

# === 3. Harmonize genes and combine matrices ===
common_genes <- intersect(rownames(counts1), rownames(counts2))
counts1 <- counts1[common_genes, , drop = FALSE]
counts2 <- counts2[common_genes, , drop = FALSE]

combined_counts <- cbind(counts2, counts1)
combined_coldata <- rbind(coldata2, coldata1)
# ensure ordering of coldata matches columns
combined_coldata <- combined_coldata[colnames(combined_counts), , drop = FALSE]

# Make sure types are correct
combined_coldata$batch <- factor(combined_coldata$batch)
combined_coldata$bio <- factor(combined_coldata$bio)
combined_coldata$treatment <- factor(combined_coldata$treatment)

# === 4. DESeq2 object with batch in design (for differential testing) ===
# Use bio (genotype) as biological factor. We do NOT force LiCl into the main contrast but keep it in colData.
dds <- DESeqDataSetFromMatrix(countData = round(as.matrix(combined_counts)),   # ensure integers
                              colData = combined_coldata,
                              design = ~ batch + bio)
dds <- DESeq(dds)   # run DESeq for DE testing (batch adjusted in the model)

# === 5. VST for visualization (preserve design by blind = FALSE) ===
vsd <- vst(dds, blind = FALSE)
mat_vst <- assay(vsd)   # genes x samples (log2-like)

# === 6. Batch correction for PCA/plots: removeBatchEffect while preserving bio + treatment ===
design_matrix <- model.matrix(~ bio + treatment, data = as.data.frame(colData(vsd)))
mat_batch_corrected <- removeBatchEffect(mat_vst, batch = combined_coldata$batch, design = design_matrix)

# Optional: if you prefer ComBat use:
# library(sva)
# modcombat <- model.matrix(~ bio + treatment, data = combined_coldata)
# mat_batch_corrected <- ComBat(dat = as.matrix(mat_vst), batch = combined_coldata$batch, mod = modcombat, par.prior = TRUE, prior.plots = FALSE)

# === 7. Filter low-variance genes before PCA (same approach as you had) ===
variances <- apply(mat_batch_corrected, 1, var)
rlog_filtered <- mat_batch_corrected[variances > 1e-6, , drop = FALSE]

# === 8. PCA: before and after correction for QC ===
pca_before <- prcomp(t(mat_vst[row.names(rlog_filtered), ]), scale. = TRUE)
pca_after  <- prcomp(t(rlog_filtered), scale. = TRUE)

explained_var_before <- round((pca_before$sdev^2 / sum(pca_before$sdev^2)) * 100, 2)
explained_var_after  <- round((pca_after$sdev^2 / sum(pca_after$sdev^2)) * 100, 2)

pca_df_before <- as.data.frame(pca_before$x)
pca_df_before$Sample <- rownames(pca_df_before)
pca_df_before$bio <- combined_coldata$bio
pca_df_before$treatment <- combined_coldata$treatment
pca_df_before$batch <- combined_coldata$batch

pca_df_after <- as.data.frame(pca_after$x)
pca_df_after$Sample <- rownames(pca_df_after)
pca_df_after$bio <- combined_coldata$bio
pca_df_after$treatment <- combined_coldata$treatment
pca_df_after$batch <- combined_coldata$batch

# === 9. Create human-friendly sample labels for plotting ===
# For Batch1 (main experiment) keep: WT1, WT2, WT3, KO1, KO2, KO3
label_map_batch1 <- setNames(colnames(counts2),
                             c("WT1","WT2","WT3","KO1","KO2","KO3"))
# For Batch2 (rescue): create KO2-1, KO2-2, KO2-3 and KO2 LiCl-1,... according to treatment
batch2_samples <- rownames(combined_coldata)[combined_coldata$batch == "Batch2"]
batch2_treatment <- combined_coldata[batch2_samples, "treatment"]
untreated_idx <- which(batch2_treatment == "none")
licl_idx <- which(batch2_treatment == "LiCl")

labels_batch2 <- character(length(batch2_samples))
labels_batch2[untreated_idx] <- paste0("KO2-", seq_along(untreated_idx))
labels_batch2[licl_idx] <- paste0("KO2 LiCl-", seq_along(licl_idx))
names(labels_batch2) <- batch2_samples

# Combine maps
sample_label_map <- c(label_map_batch1, labels_batch2)
# If any sample names aren't in the map (rare), default to the raw column name
all_samples <- colnames(combined_counts)
sample_labels <- ifelse(all_samples %in% names(sample_label_map),
                        sample_label_map[all_samples],
                        all_samples)

# Add friendly label to PCA dfs (respect rownames)
pca_df_before$Label <- sample_labels[match(pca_df_before$Sample, all_samples)]
pca_df_after$Label  <- sample_labels[match(pca_df_after$Sample, all_samples)]

# For plotting group (color): we want WT, KO1, KO2, KO3, and KO2+LiCl
pca_df_before$Group <- ifelse(pca_df_before$treatment == "LiCl", "KO2+LiCl", as.character(pca_df_before$bio))
pca_df_after$Group  <- ifelse(pca_df_after$treatment == "LiCl", "KO2+LiCl", as.character(pca_df_after$bio))

# === 10. Colors ===
group_colors <- c(
  "WT" = "#4DAF4A",
  "KO1" = "#E41A1C",
  "KO2" = "#984EA3",
  "KO3" = "#FF7F00",
  "KO2+LiCl" = "#377EB8"
)


# --- Rename specific samples for plotting (insert just before plotting) ---
# Map: original column name -> plotting label
rename_map <- c(
  "KO_75" = "KO2 batch1",
  "WT_1"  = "WT1",
  "WT_2"  = "WT2",
  "WT_3"  = "WT3",
  "KO_9"  = "KO1",
  "KO_46" = "KO3"
)

# all_samples should be colnames(combined_counts) from earlier
# Update the sample_label_map (if present) and sample_labels vector used for plotting
if (exists("sample_label_map")) {
  for (orig in names(rename_map)) {
    if (orig %in% names(sample_label_map)) sample_label_map[orig] <- rename_map[orig]
  }
}
if (exists("all_samples") && exists("sample_labels")) {
  for (orig in names(rename_map)) {
    idx <- which(all_samples == orig)
    if (length(idx) == 1) sample_labels[idx] <- rename_map[orig]
  }
}

# Re-assign the friendly labels in the PCA data frames (fall back to the raw sample name if no previous label)
pca_df_before$Label <- sample_labels[match(pca_df_before$Sample, all_samples)]
pca_df_after$Label  <- sample_labels[match(pca_df_after$Sample, all_samples)]

# As a safety: if any Sample still didn't get a friendly label, replace with the raw sample name
pca_df_before$Label[is.na(pca_df_before$Label) | pca_df_before$Label == ""] <- pca_df_before$Sample[is.na(pca_df_before$Label) | pca_df_before$Label == ""]
pca_df_after$Label[is.na(pca_df_after$Label)  | pca_df_after$Label == ""]  <- pca_df_after$Sample[is.na(pca_df_after$Label)  | pca_df_after$Label == ""]

# (Optional) Print mapping applied for quick check
message("Applied plotting renames:\n", paste(paste(names(rename_map), "->", rename_map), collapse = "\n"))

# === 11. PCA plots: before and after ===
p_before <- ggplot(pca_df_before, aes(x = PC1, y = PC2, color = Group, label = Label)) +
  geom_point(size = 4, alpha = 0.9) +
  geom_text_repel(size = 4, show.legend = FALSE) +
  scale_color_manual(values = group_colors) +
  labs(title = paste0("PCA BEFORE batch correction (PC1 ", explained_var_before[1], "%; PC2 ", explained_var_before[2], "%)"),
       x = paste0("PC1 (", explained_var_before[1], "%)"),
       y = paste0("PC2 (", explained_var_before[2], "%)")) +
  theme_minimal(base_size = 14) + theme(legend.position = "top", plot.title = element_text(hjust = 0.5))

p_after <- ggplot(pca_df_after, aes(x = PC1, y = PC2, color = Group, label = Label)) +
  geom_point(size = 4, alpha = 0.9) +
  geom_point(data = pca_df_after[pca_df_after$Group == "KO2", ], size = 6, shape = 21, fill = "#984EA3", color = "#984EA3") +
  geom_text_repel(size = 4, show.legend = FALSE) +
  scale_color_manual(values = group_colors) +
  labs(title = paste0("PCA AFTER removeBatchEffect (PC1 ", explained_var_after[1], "%; PC2 ", explained_var_after[2], "%)"),
       x = paste0("PC1 (", explained_var_after[1], "%)"),
       y = paste0("PC2 (", explained_var_after[2], "%)")) +
  theme_minimal(base_size = 14) + theme(legend.position = "top", plot.title = element_text(hjust = 0.5))

# Print both plots
print(p_before)
print(p_after)

