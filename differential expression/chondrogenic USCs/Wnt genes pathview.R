#KEGG pathway diagrams top20 pathways

# Step 1: Load necessary libraries
library(biomaRt)
library(clusterProfiler)
library(pathview)
library(org.Hs.eg.db)
library(readxl)

# Step 2: Read in your DEG data from Excel
deg_file <- "W:/humangenetik/data/AG Thiel/Mitarbeiter/Alexander/LAMA5/CRISPR/USC (LAMA5 KO)/KOs/LAMA5 KO C Nr.2 (FACS) - Erfolgreich/RNA Seq/ZelltypenEinzeln WT vs. KO/differentialExpression_Thiel_Chond.xlsx"
deg_table <- read_excel(deg_file)

# Step 3: Get Entrez IDs for the DEG table genes
# Mapping ENSEMBL IDs to Entrez IDs
ensembl_ids <- deg_table$ENSEMBL
entrez_ids <- mapIds(org.Hs.eg.db, keys = ensembl_ids, column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first")

# Add Entrez IDs to the DEG table
deg_table$entrez_id <- entrez_ids

# Step 4: Perform KEGG pathway enrichment analysis on all DEGs
deg_entrez_ids <- na.omit(deg_table$entrez_id)  # Ensure no NA values are passed
pathway_enrichment <- enrichKEGG(gene = deg_entrez_ids, organism = "hsa")

# Step 1: Sort pathways by FoldEnrichment to select top 20 pathways
top_pathways <- pathway_enrichment@result[order(-as.numeric(pathway_enrichment@result$FoldEnrichment)), ]
top_20_pathways <- head(top_pathways, 20)  # Select top 20 pathways based on FoldEnrichment

# Step 2: Extract the list of gene IDs for the top 20 pathways
top_20_pathways_genes <- list()
for (i in 1:nrow(top_20_pathways)) {
  pathway_id <- top_20_pathways$ID[i]
  genes <- top_20_pathways$geneID[i]
  top_20_pathways_genes[[pathway_id]] <- strsplit(genes, "/")[[1]]
}

# Step 3: Map the gene IDs from the DEG table
# We already have the DEG table with Entrez IDs, so no need to do mapping again.
# Ensure we only include genes present in the DEG table
pathway_genes_deg <- data.frame()
for (pathway_id in names(top_20_pathways_genes)) {
  pathway_genes <- top_20_pathways_genes[[pathway_id]]
  deg_for_pathway <- deg_table[deg_table$entrez_id %in% pathway_genes, c("entrez_id", "log2FoldChange")]
  deg_for_pathway$pathway <- pathway_id
  pathway_genes_deg <- rbind(pathway_genes_deg, deg_for_pathway)
}

# Step 4: Revert log2FoldChange (KO vs WT) for consistency
pathway_genes_deg$log2FoldChange <- -pathway_genes_deg$log2FoldChange

# Step 5: Manually scale log2FoldChange to range [-5, 5]
manual_scale <- function(x) {
  x <- pmin(pmax(x, -5), 5)
  return(x)
}

# Step 6: Annotate the log2FoldChange in pathway diagrams and save to default location (Documents)
for (pathway in names(top_20_pathways_genes)) {
  
  # Extract pathway data for the current pathway
  pathway_data <- pathway_genes_deg[pathway_genes_deg$pathway == pathway, ]
  pathway_entrez_ids <- pathway_data$entrez_id
  pathway_log2fc <- pathway_data$log2FoldChange
  
  # Filter out genes with missing log2FoldChange values
  pathway_data_filtered <- pathway_data[!is.na(pathway_data$log2FoldChange), ]
  pathway_entrez_ids_filtered <- pathway_data_filtered$entrez_id
  pathway_log2fc_filtered <- pathway_data_filtered$log2FoldChange
  
  # Check if there are any valid genes left after filtering
  if (length(pathway_entrez_ids_filtered) == 0) {
    print(paste("Skipping pathway", pathway, "because all log2FoldChange values are missing"))
    next  # Skip this pathway if all genes have missing data
  }
  
  # Apply the manual scaling to the log2FoldChange values
  pathway_log2fc_scaled <- manual_scale(pathway_log2fc_filtered)
  
  # Create a named vector for log2FC to annotate the pathway diagram
  log2fc_vector <- setNames(pathway_log2fc_scaled, as.character(pathway_entrez_ids_filtered))
  
  # Define the file path for saving the diagram in default directory (Documents)
  output_file <- paste0("pathway_", pathway, ".png")  # Save directly to current working directory
  
  # Print output_file to debug the file path
  print(paste("Saving to:", output_file))
  
  # Ensure pathview function saves the file correctly
  png(output_file, width = 800, height = 600)  # Set desired resolution
  tryCatch({
    pv <- pathview(
      gene.data = log2fc_vector,
      pathway.id = pathway,
      species = "hsa",
      cex.lab = 0.5,
      cex.legend = 0.5,
      low = "blue",
      high = "red",
      na.col = "white",
      legend.position = "bottomright",
      legend.title = "log2FC",
      legend.cex = 0.8,
      show.cnames = FALSE
    )
  }, error = function(e) {
    print(paste("Error while processing pathway", pathway, ":", e$message))
  })
  
  dev.off()
  
  print(paste("Successfully saved:", output_file))
}

