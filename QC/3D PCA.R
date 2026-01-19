
# Load required libraries
library(DESeq2)
library(ggplot2)
library(plotly)
library(readxl)
library(htmlwidgets)
library(alphashape3d) # For alpha shapes

# Load and process data
filepath <- "...Counts_ThielRNA.xlsx"
count_data <- read_excel(filepath, sheet = 1)
counts <- as.data.frame(count_data[,-1])
rownames(counts) <- count_data$ID
condition <- factor(c(rep("WT", 3), rep("KO", 3), rep("CHONDR_WT", 3), rep("CHONDR_KO", 3),
                      rep("OsteoLAMA5_WT", 3), rep("OsteoLAMA5_KO", 3)))
coldata <- data.frame(row.names = colnames(counts), condition)
dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ condition)
normalized_counts <- assay(rlog(dds))

# Filter low variance and constant rows/columns
variance_threshold <- 1e-6
variances <- apply(normalized_counts, 2, var)
low_variance_cols <- variances < variance_threshold
normalized_counts_filtered <- normalized_counts[, !low_variance_cols]
constant_cols <- apply(normalized_counts_filtered, 2, function(x) var(x) == 0)
normalized_counts_filtered <- normalized_counts_filtered[, !constant_cols]
constant_rows <- apply(normalized_counts_filtered, 1, function(x) var(x) == 0)
normalized_counts_filtered <- normalized_counts_filtered[!constant_rows, ]

# Perform PCA
pca_res <- prcomp(t(normalized_counts_filtered), scale. = TRUE)
pca_data <- as.data.frame(pca_res$x)

# Rename samples for clarity
new_names <- c("WT1", "WT2", "WT3", "KO1", "KO2", "KO3", 
               "chWT1", "chWT2", "chWT3", "chKO1", "chKO2", "chKO3",
               "oWT1", "oWT2", "oWT3", "oKO1", "oKO2", "oKO3")
colnames(normalized_counts_filtered) <- new_names
filtered_samples <- colnames(normalized_counts_filtered)
pca_data$Sample <- filtered_samples

# Add group information and colors
group_mapping <- c(rep("WT", 3), rep("KO", 3), rep("chWT", 3), rep("chKO", 3), 
                   rep("oWT", 3), rep("oKO", 3))
pca_data$Group <- factor(group_mapping, levels = unique(group_mapping))

# Define high-contrast colors manually
group_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#A65628") 

# Initialize 3D PCA plot
fig <- plot_ly()

# Add points for each group
for (grp in levels(pca_data$Group)) {
  group_points <- subset(pca_data, Group == grp)
  fig <- fig %>%
    add_markers(
      data = group_points,
      x = ~PC1, y = ~PC2, z = ~PC3,
      color = I(group_colors[which(levels(pca_data$Group) == grp)]),
      name = grp,
      marker = list(size = 6, opacity = 0.8)
    )
}

# Add triangular meshes for each group (assuming 3 samples per group)
for (grp in levels(pca_data$Group)) {
  group_points <- subset(pca_data, Group == grp)[, c("PC1", "PC2", "PC3")]
  
  if (nrow(group_points) == 3) {
    fig <- fig %>%
      add_mesh(
        x = group_points$PC1,
        y = group_points$PC2,
        z = group_points$PC3,
        opacity = 0.2,
        color = I(group_colors[which(levels(pca_data$Group) == grp)]),
        name = paste(grp, "Triangle")
      )
  }
}

# Customize layout
fig <- fig %>%
  layout(
    title = list(text = "3D PCA Plot with Group Shapes", font = list(size = 20)),
    scene = list(
      xaxis = list(title = list(text = "PC1", font = list(size = 16))),
      yaxis = list(title = list(text = "PC2", font = list(size = 16))),
      zaxis = list(title = list(text = "PC3", font = list(size = 16)))
    ),
    legend = list(title = list(text = "Groups", font = list(size = 16))),
    font = list(size = 14)
  )

# Export the plot
saveWidget(fig, "pca_plot_with_alpha_shapes.html")
