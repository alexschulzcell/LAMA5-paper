# Load libraries
library(DESeq2)
library(readxl)
library(writexl)
library(org.Hs.eg.db)

# === 1. Load counts ===
file_path <- "...ThielRNA.tab"

count_data <- read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Clean sample names
colnames(count_data) <- gsub("^ThielRNA_|\\.bam$", "", colnames(count_data))

# Set gene IDs as rownames
rownames(count_data) <- count_data$Geneid
counts <- count_data[, -1]

# === 2. Define conditions ===
sample_conditions <- factor(c(rep("LiCl", 3), rep("KO", 3)))
coldata <- data.frame(condition = sample_conditions, row.names = colnames(counts))

# === 3. Create DESeq2 dataset ===
dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ condition)

# Run DESeq2
dds <- DESeq(dds)

# === 4. Get results ===
res <- results(dds)

# Filter out NA padj values
res_filtered <- res[!is.na(res$padj), ]

# Filter significant results (padj < 0.05)
res_sig <- res_filtered[res_filtered$padj < 0.05, ]

# Convert to data frame and add gene IDs
res_filtered_df <- as.data.frame(res_filtered)
res_filtered_df$gene <- rownames(res_filtered_df)

# Map Ensembl IDs to gene symbols
gene_symbols <- mapIds(org.Hs.eg.db, 
                       keys = res_filtered_df$gene, 
                       column = "SYMBOL", 
                       keytype = "ENSEMBL", 
                       multiVals = "first")
res_filtered_df$gene_symbol <- gene_symbols

# === 5. Write results to Excel ===
write_xlsx(list("DESeq2 Results" = res_filtered_df), path = "DESeq2_KO_vs_LiCl_results.xlsx")
