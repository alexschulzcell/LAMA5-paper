# Load required libraries
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(htmlwidgets)

# === 1. Load data ===
filepath <- "...ThielRNA.tab"

# Read tab-separated file
count_data <- read.table(filepath, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# First column = gene IDs
rownames(count_data) <- count_data$Geneid
counts <- count_data[, -1]

# === 2. Define conditions ===
condition <- factor(c(rep("KO", 3), rep("LiCl", 3)))
coldata <- data.frame(row.names = colnames(counts), condition)

# === 3. Create DESeq2 dataset and normalize ===
dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ condition)
dds <- DESeq(dds)
normalized_counts <- assay(rlog(dds, blind = TRUE))

# === 4. Filter out low-variance and constant genes ===
variances <- apply(normalized_counts, 1, var)
normalized_counts_filtered <- normalized_counts[variances > 1e-6, ]

# === 5. PCA ===
pca_res <- prcomp(t(normalized_counts_filtered), scale. = TRUE)
pca_data <- as.data.frame(pca_res$x)
pca_data$Sample <- rownames(pca_data)
pca_data$Group <- condition

# Calculate explained variance for axis labels
explained_var <- round((pca_res$sdev^2 / sum(pca_res$sdev^2)) * 100, 2)
pca_data$Sample <- gsub("^ThielRNA_|\\.bam$", "", pca_data$Sample)



# === 6. Define colors ===
group_colors <- c("KO" = "#E41A1C", "LiCl" = "#377EB8")

# === 7. 2D PCA plot ===
p <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Group, label = Sample)) +
  geom_point(size = 5, alpha = 0.8) +
  geom_text_repel(size = 5, show.legend = FALSE) +
  scale_color_manual(values = group_colors) +
  labs(
    title = "PCA Plot (KO vs LiCl)",
    x = paste0("PC1 (", explained_var[1], "%)"),
    y = paste0("PC2 (", explained_var[2], "%)"),
    color = "Group"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "top",
    legend.title = element_text(face = "bold")
  )

print (p)
