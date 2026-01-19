# ===============================
# Human Fetal Growth Plate scRNA
# LAMA5 – PITX1 – FLI1 Analysis
# ===============================

library(Seurat)
library(SeuratDisk)
library(dplyr)
library(ggplot2)
library(patchwork)

# -------------------------------
# 1. Load and Merge Loom Files
# -------------------------------

loom_files <- c(
  "...5478STDY7980349.loom",
  "...5478STDY7935101.loom",
  "...5478STDY7935102.loom",
  "...5386STDY7537944.loom",
  "...5478STDY7980348.loom"
)

seurat_list <- lapply(seq_along(loom_files), function(i) {
  seu <- LoadLoom(loom_files[i])
  seu$sample_id <- paste0("HL_", i)
  seu
})

seu <- merge(
  seurat_list[[1]],
  y = seurat_list[-1],
  add.cell.ids = paste0("HL", seq_along(seurat_list))
)

# -------------------------------
# 2. Standard Preprocessing
# -------------------------------

seu <- NormalizeData(seu, verbose = FALSE)
seu <- FindVariableFeatures(seu, nfeatures = 3000, verbose = FALSE)
seu <- ScaleData(seu, verbose = FALSE)

seu <- RunPCA(seu, npcs = 30, verbose = FALSE)
seu <- RunUMAP(seu, dims = 1:30, verbose = FALSE)

seu <- FindNeighbors(seu, dims = 1:30, verbose = FALSE)
seu <- FindClusters(seu, resolution = 0.5, verbose = FALSE)

DimPlot(seu, label = TRUE) + ggtitle("UMAP – Unannotated Clusters")

# -------------------------------
# STEP 3 – Automatic Cell-Type Annotation (ROBUST)
# -------------------------------

# Define marker sets
marker.list <- list(
  Chondrocytes = c("SOX9", "COL2A1", "ACAN"),
  Perichondrium = c("PRRX1", "COL1A1", "TWIST1"),
  Endothelium = c("PECAM1", "KDR", "EMCN")
)

# Add module scores
seu <- AddModuleScore(
  seu,
  features = marker.list,
  name = names(marker.list)
)

# Identify newly created columns
score.cols <- grep(
  paste0("^", names(marker.list), collapse = "|"),
  colnames(seu@meta.data),
  value = TRUE
)

# Sanity check
print(score.cols)

# Create mapping (ordered by marker.list)
score.map <- setNames(score.cols, names(marker.list))

# Summarise per cluster
cluster_scores <- seu@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(
    Chondrocytes = mean(.data[[score.map["Chondrocytes"]]]),
    Perichondrium = mean(.data[[score.map["Perichondrium"]]]),
    Endothelium = mean(.data[[score.map["Endothelium"]]])
  )

# Assign dominant identity
cluster_scores$celltype <- apply(
  cluster_scores[, c("Chondrocytes", "Perichondrium", "Endothelium")],
  1,
  function(x) names(x)[which.max(x)]
)

# Apply identities
new.cluster.ids <- setNames(
  cluster_scores$celltype,
  cluster_scores$seurat_clusters
)

seu <- RenameIdents(seu, new.cluster.ids)
seu$celltype <- Idents(seu)

# UMAP with auto-annotated labels
DimPlot(
  seu,
  group.by = "celltype",
  label = TRUE,
  repel = TRUE
) +
  ggtitle("Automatically Annotated Cell Types")

# Marker confirmation
DotPlot(
  seu,
  features = c("SOX9", "COL2A1", "ACAN", "PRRX1", "COL1A1", "TWIST1", "PECAM1", "KDR", "EMCN"),
  group.by = "celltype"
) +
  RotatedAxis()

DimPlot(
  seu,
  group.by = "celltype",
  label = TRUE,
  repel = TRUE
) +
  ggtitle("Panel A – Cell-Type Annotation of Human Fetal Growth Plate") +
  theme(legend.position = "right")


markers.to.plot <- c(
  "SOX9", "COL2A1", "ACAN", "PRRX1", "COL1A1", "TWIST1",
  "PECAM1", "KDR", "EMCN", "FLI1",
  "PITX1", "WNT7A", "TFAP2A", "GRHL2", "LAMA5"
)

dp <- DotPlot(
  seu,
  features = markers.to.plot,
  group.by = "celltype"
)

dp_data <- dp$data

ggplot(dp_data, aes(x = id, y = features.plot)) +
  geom_point(
    aes(size = pct.exp, fill = avg.exp.scaled),
    shape  = 21,
    color  = "black",
    stroke = 0.3
  ) +
  scale_fill_gradientn(
    colors = c("#66A3FF", "white", "#F66666"),
    name   = "Expression"
  ) +
  scale_size(range = c(1, 8), name = "Pct. Expressing") +
  labs(
    title = "Dotplot of Network Genes",
    x     = NULL,
    y     = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title       = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.text.x      = element_text(size = 14, angle = 45, hjust = 1),
    axis.text.y      = element_text(size = 14, face = "italic"),
    legend.title     = element_text(size = 14),
    legend.text      = element_text(size = 12),
    panel.grid.major = element_line(color = "grey90")
  )

markers.to.plot <- c(
  # 1) Sender: ECM + ligand
  "LAMA5", "WNT7A",
  # 2) Boundary
  "TFAP2A", "GRHL2",
  # 3) Receiver: transcriptional targets
  "PITX1", "FLI1",
  # 4) Identity – chondrocytes first (receiver), then perichondrium (sender base)
  "SOX9", "COL2A1", "ACAN",
  "PRRX1", "COL1A1", "TWIST1",
  # 5) Vascular reference
  "PECAM1", "KDR", "EMCN"
)

dp <- DotPlot(
  seu,
  features = markers.to.plot,
  group.by = "celltype"
)

dp_data <- dp$data

# Order cell types along your axis: Sender -> Receiver -> Reference
dp_data$id <- factor(
  dp_data$id,
  levels = c("Perichondrium", "Chondrocytes", "Endothelium")
)

# Keep gene order as defined in markers.to.plot
dp_data$features.plot <- factor(
  dp_data$features.plot,
  levels = rev(markers.to.plot)  # rev() so top row = first element in vector
)

ggplot(dp_data, aes(x = id, y = features.plot)) +
  geom_point(
    aes(size = pct.exp, fill = avg.exp.scaled),
    shape  = 21,
    color  = "black",
    stroke = 0.3
  ) +
  scale_fill_gradientn(
    colors = c("#66A3FF", "white", "#F66666"),
    name   = "Expression"
  ) +
  scale_size(range = c(1, 8), name = "Pct. Expressing") +
  labs(
    title = "Perichondrial Sender – Chondrocyte Receiver – Vascular Reference",
    x     = NULL,
    y     = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title       = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.text.x      = element_text(size = 18, angle = 45, hjust = 1),
    axis.text.y      = element_text(size = 18, face = "italic"),
    legend.title     = element_text(size = 16),
    legend.text      = element_text(size = 16),
    panel.grid.major = element_line(color = "grey90")
  )

# ==============================================================================
# FIGURE 4B: The Complete Niche Map (LAMA5 + Signal + Target)
# ==============================================================================

# 1. Define the complete gene list
genes_complete <- c("LAMA5", "WNT7A",  # The Niche / Matrix / Signal
                    "TFAP2A", "GRHL2", # Perichondrial Identity
                    "FLI1", "PITX1",   # The Chondrocytescyte Targets
                    "SOX9", "COL1A1")  # Broad Identity

# Fetch data
expr_data <- FetchData(seu, vars = genes_complete)

# --------------------------------------------------
# 0. Prepare expression matrix for the key genes
#    (rows = cells, cols = genes)
# --------------------------------------------------

genes_use <- c("LAMA5", "WNT7A", "TFAP2A", "GRHL2",
               "PITX1", "FLI1")

expr_data <- t(as.matrix(
  GetAssayData(seu, assay = "RNA", slot = "data")[genes_use, ]
))

# Sanity: same order
stopifnot(nrow(expr_data) == nrow(seu@meta.data))

celltype <- as.character(seu$celltype)

# --------------------------------------------------
# 1. Define the Hierarchy (Niche -> Signal -> Target)
#    with background = existing celltype
# --------------------------------------------------

seu$Niche_Architecture <- vapply(
  X = seq_len(nrow(expr_data)),
  FUN = function(i) {
    x <- expr_data[i, ]
    
    # --- GROUP A: THE ARCHITECTS (Perichondrial/Niche) ---
    
    # 1. The "Holy Grail": Matrix + Ligand
    if (x["LAMA5"] > 0 && x["WNT7A"] > 0) {
      return("Niche Hub (LAMA5+/WNT7A+)")
    }
    
    # 2. The Mechanics: Matrix-only layer
    if (x["LAMA5"] > 0) {
      return("Matrix Architect (LAMA5+)")
    }
    
    # 3. The Signal: Ligand-only (non‑LAMA5)
    if (x["WNT7A"] > 0) {
      return("Niche Signal (WNT7A+)")
    }
    
    # 4. The Boundary: Perichondrial module
    if (x["TFAP2A"] > 0 || x["GRHL2"] > 0) {
      return("Boundary (TFAP2A+/GRHL2+)")
    }
    
    # --- GROUP B: THE TARGETS (Chondrocytes) ---
    
    # 5. Co‑regulated targets
    if (x["PITX1"] > 0 && x["FLI1"] > 0) {
      return("Target: Dual Positive (PITX1+/FLI1+)")
    }
    if (x["PITX1"] > 0) {
      return("Target: PITX1+")
    }
    if (x["FLI1"] > 0) {
      return("Target: FLI1+")
    }
    
    # --- GROUP C: BACKGROUND ---
    # For all remaining cells, just use the coarse celltype
    # (e.g. "Chondrocytes", "Perichondrium", "Endothelium")
    return(celltype[i])
  },
  FUN.VALUE = character(1)
)



plot_order <- c(
  "Niche Hub (LAMA5+/WNT7A+)",
  "Niche Signal (WNT7A+)",
  "Matrix Architect (LAMA5+)",
  "Boundary (TFAP2A+/GRHL2+)",
  "Target: Dual Positive (PITX1+/FLI1+)",
  "Target: PITX1+",
  "Target: FLI1+",
  "Chondrocytes",
  "Perichondrium",
  "Endothelium"
)

niche_colors <- c(
  "Niche Hub (LAMA5+/WNT7A+)"            = "green",
  "Niche Signal (WNT7A+)"                = "#4B0092",
  "Matrix Architect (LAMA5+)"            = "#332288",
  "Boundary (TFAP2A+/GRHL2+)"            = "#88CCEE",
  "Target: Dual Positive (PITX1+/FLI1+)" = "#D55E00",
  "Target: PITX1+"                       = "#E69F00",
  "Target: FLI1+"                        = "#CC79A7",
  "Chondrocytes"                         = "#F0E442",
  "Perichondrium"                        = "#DDDDDD",
  "Endothelium"                          = "grey70",
  "Other"                                = "grey98"
)




p <- DimPlot(
  seu,
  group.by = "Niche_Architecture",
  cols     = niche_colors,
  order    = plot_order,
  pt.size  = 1.2
) +
  ggtitle("Human Signaling Niche") +
  theme_void() +
  theme(
    plot.title      = element_text(size = 18, face = "bold", hjust = 0.5),
    legend.position = "right",
    legend.title    = element_blank(),
    legend.text     = element_text(size = 18)
  )

umap_mat <- Embeddings(seu, "umap")
xr <- range(umap_mat[,1])
yr <- range(umap_mat[,2])

dx <- diff(xr) * 0.20
dy <- diff(yr) * 0.20

p +
  geom_segment(
    aes(x = xr[1], y = yr[1], xend = xr[1] + dx, yend = yr[1]),
    arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
    linewidth = 0.8,
    inherit.aes = FALSE
  ) +
  geom_segment(
    aes(x = xr[1], y = yr[1], xend = xr[1], yend = yr[1] + dy),
    arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
    linewidth = 0.8,
    inherit.aes = FALSE
  ) +
  annotate(
    "text",
    x = xr[1],
    y = yr[1] - 0.02 * diff(yr),
    label = "UMAP1",
    hjust = 0,
    vjust = 1,
    size = 6
  ) +
  annotate(
    "text",
    x = xr[1] - 0.02 * diff(xr),
    y = yr[1],
    label = "UMAP2",
    hjust = 0,
    vjust = 0,
    angle = 90,
    size = 6
  )


celltype_colors <- c(
  "Chondrocytes"  = "#F0E442",
  "Perichondrium" = "#DDDDDD",
  "Endothelium"   = "grey70"
)


# -------------------------------
# Cell type UMAP with same design
# -------------------------------

ct_p <- DimPlot(
  seu,
  group.by   = "celltype",
  cols       = celltype_colors,
  label      = FALSE
) +
  ggtitle("Cell Type Annotation") +
  theme_void() +
  theme(
    plot.title      = element_text(size = 18, face = "bold", hjust = 0.5),
    legend.position = "right",
    legend.title    = element_blank(),
    legend.text     = element_text(size = 18)
  )

# Reuse same UMAP coordinate ranges for the arrows
umap_mat <- Embeddings(seu, "umap")
xr <- range(umap_mat[, 1])
yr <- range(umap_mat[, 2])

dx <- diff(xr) * 0.20
dy <- diff(yr) * 0.20

ct_p_with_axes <- ct_p +
  geom_segment(
    aes(x = xr[1], y = yr[1], xend = xr[1] + dx, yend = yr[1]),
    arrow       = arrow(length = unit(0.25, "cm"), type = "closed"),
    linewidth   = 0.9,
    inherit.aes = FALSE
  ) +
  geom_segment(
    aes(x = xr[1], y = yr[1], xend = xr[1], yend = yr[1] + dy),
    arrow       = arrow(length = unit(0.25, "cm"), type = "closed"),
    linewidth   = 0.9,
    inherit.aes = FALSE
  ) +
  annotate(
    "text",
    x     = xr[1],
    y     = yr[1] - 0.02 * diff(yr),
    label = "UMAP1",
    hjust = 0,
    vjust = 1,
    size  = 6
  ) +
  annotate(
    "text",
    x     = xr[1] - 0.02 * diff(xr),
    y     = yr[1],
    label = "UMAP2",
    hjust = 0,
    vjust = 0,
    angle = 90,
    size  = 6
  )

ct_p_with_axes

