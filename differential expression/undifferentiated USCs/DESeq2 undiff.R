DeSeq2 USCs undifferenziert WT vs. KO
library(DESeq2)
library(readxl)
library(writexl) 
library(org.Hs.eg.db)

# Load data
file_path <- "...Counts_ThielRNA.xlsx"

# Load the raw counts data from sheet 1
raw_counts <- read_excel(file_path, sheet = 1)

# Convert raw_counts to a data.frame (this avoids the tibble warning)
raw_counts <- as.data.frame(raw_counts)

# Set the gene IDs (first column) as rownames of the raw counts matrix
rownames(raw_counts) <- raw_counts$ID  # Assuming 'ID' is the first column with gene IDs
raw_counts <- raw_counts[, -1]  # Remove the 'ID' column from the matrix

# Subset the raw counts to only include the relevant samples (excluding KO samples)
raw_counts_subset <- raw_counts[, c("WT_1", "WT_2", "WT_3", "KO_9", "KO_75", "KO_46")]

# Create the sample metadata (condition for DESeq2)
sample_conditions <- factor(c("WT", "WT", "WT", "KO", "KO", "KO"))
coldata <- data.frame(condition = sample_conditions)

# Create DESeq2 dataset using raw counts with gene IDs as rownames
dds <- DESeqDataSetFromMatrix(countData = raw_counts_subset, colData = coldata, design = ~ condition)

# Convert counts to integer mode (DESeq2 requires counts to be integers)
dds <- DESeq(dds)

# Get results (differential expression)
res <- results(dds)

# Filter significant results (padj < 0.05)
res_sig <- res[which(res$padj < 0.05),]

# Filter out NA values from DESeq2 results
res_filtered <- res[!is.na(res$padj), ]

# Convert DESeq2 results to a data frame and add gene column (row names)
res_filtered_df <- as.data.frame(res_filtered)
res_filtered_df$gene <- rownames(res_filtered)

# Map gene IDs to gene symbols using org.Hs.eg.db
gene_symbols <- mapIds(org.Hs.eg.db, keys = res_filtered_df$gene, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")

# Add gene symbols as a new column to the DESeq2 results
res_filtered_df$gene_symbol <- gene_symbols

# Write to Excel
write_xlsx(list("DESeq2 Results" = res_filtered_df), path = "DESeq2_results.xlsx")
