library(data.table)
library(Seurat)
library(dplyr)

raw_dir <- "...GSE185940_RAW"
files   <- list.files(raw_dir, pattern = "\\.txt\\.gz$", full.names = TRUE)

length(files)  # should be 103
read_gsm_matrix <- function(f) {
  message("Reading: ", basename(f))
  dt <- fread(f)
  
  # First column is genes (currently called V1)
  gene_col <- colnames(dt)[1]  # "V1"
  genes    <- dt[[gene_col]]
  
  mat <- as.matrix(dt[, -1, with = FALSE])
  rownames(mat) <- genes
  
  # Make cell names unique by prefixing with the GSM ID
  sample_id <- tools::file_path_sans_ext(basename(f))  # GSM5627543_AB4005
  colnames(mat) <- paste0(sample_id, "_", colnames(mat))
  
  mat
}
mat_list <- lapply(files, read_gsm_matrix)

# Quick sanity check: number of genes in first few matrices
sapply(mat_list[1:3], nrow)
# Compare the first file's gene list to a few others
g1 <- rownames(mat_list[[1]])

all_same <- all(sapply(mat_list, function(m) identical(rownames(m), g1)))
all_same
combined_mat <- do.call(cbind, mat_list)
dim(combined_mat)  # genes × cells


seu_m <- CreateSeuratObject(
  counts       = combined_mat,
  project      = "GSE185940_mouse_limb",
  min.cells    = 3,
  min.features = 200
)

seu_m


seu_m <- NormalizeData(seu_m)
seu_m <- FindVariableFeatures(seu_m, selection.method = "vst", nfeatures = 2000)
seu_m <- ScaleData(seu_m)
seu_m <- RunPCA(seu_m, npcs = 30)
seu_m <- FindNeighbors(seu_m, dims = 1:20)
seu_m <- FindClusters(seu_m, resolution = 0.6)
seu_m <- RunUMAP(seu_m, dims = 1:20)
meta_path <- "...GSE185940_M_metadata_s.txt.gz"

# Open gz file as text connection
con <- gzfile(meta_path, open = "rt")

# First line is that long header string, we just read and ignore for now
header_raw <- readLines(con, n = 1)
header_raw
# "Well_ID.well_coordinates.Amp_batch_ID.Cell_barcode.Seq_batch_ID.Pool_barcode.Number_of_cells.total.genes"

# Now read the rest, splitting on *whitespace* (not tabs)
meta <- read.table(
  con,
  header           = FALSE,
  sep              = "",
  stringsAsFactors = FALSE,
  quote            = "",
  comment.char     = ""
)

close(con)

dim(meta)
head(meta)

colnames(meta) <- c(
  "Well_ID",        # 1, 2, 3, ...
  "Cell_barcode",   # WMC1072673, ...
  "Well_coord",     # A1, C1, ...
  "Amp_batch_ID",   # AB6178, ...
  "Seq_batch_ID",   # CTATTCG, ...
  "Pool_barcode",   # SB308, ...
  "Index_seq",      # CGTACCGC (second barcode)
  "Number_of_cells",
  "total_counts",
  "n_genes"
)

dim(meta)
head(meta)

# Seurat cell names you created, e.g. "GSM5627543_AB4005_WMC241697"
head(colnames(seu_m))

# extract the WMC part at the end
barcode_seurat <- sub(".*_", "", colnames(seu_m))
head(barcode_seurat)

# Make sure meta barcodes are character
meta$Cell_barcode <- as.character(meta$Cell_barcode)

# Match Seurat cells to metadata
match_idx <- match(barcode_seurat, meta$Cell_barcode)

table(is.na(match_idx))  # how many cells didn't find metadata?

# Order metadata to align with Seurat columns
meta_ordered <- meta[match_idx, ]
rownames(meta_ordered) <- colnames(seu_m)

# Attach
seu_m <- AddMetaData(seu_m, metadata = meta_ordered)




library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
# Mouse limb Seurat object from GSE185940
# (created from the combined GSM matrices)
seu_m
# -------------------------------
# STEP 1 – Automatic Cell-Type Annotation (Mouse)
# -------------------------------

marker.list.m <- list(
  Chondrocytes  = c("Sox9", "Col2a1", "Acan"),
  Perichondrium = c("Prrx1", "Col1a1", "Twist1"),
  Endothelium   = c("Pecam1", "Kdr", "Emcn")
)

seu_m <- AddModuleScore(
  seu_m,
  features = marker.list.m,
  name     = names(marker.list.m)
)

# Identify the created score columns
score.cols.m <- grep(
  paste0("^", names(marker.list.m), collapse = "|"),
  colnames(seu_m@meta.data),
  value = TRUE
)
print(score.cols.m)

score.map.m <- setNames(score.cols.m, names(marker.list.m))

cluster_scores_m <- seu_m@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(
    Chondrocytes  = mean(.data[[score.map.m["Chondrocytes"]]]),
    Perichondrium = mean(.data[[score.map.m["Perichondrium"]]]),
    Endothelium   = mean(.data[[score.map.m["Endothelium"]]])
  )

cluster_scores_m$celltype <- apply(
  cluster_scores_m[, c("Chondrocytes", "Perichondrium", "Endothelium")],
  1,
  function(x) names(x)[which.max(x)]
)

new.cluster.ids.m <- setNames(
  cluster_scores_m$celltype,
  cluster_scores_m$seurat_clusters
)

seu_m <- RenameIdents(seu_m, new.cluster.ids.m)
seu_m$celltype <- Idents(seu_m)
DimPlot(
  seu_m,
  group.by = "celltype",
  label    = TRUE,
  repel    = TRUE
) +
  ggtitle("Mouse – Marker-guided Cell Types")
DotPlot(
  seu_m,
  features = c(
    "Sox9", "Col2a1", "Acan",
    "Prrx1", "Col1a1", "Twist1",
    "Pecam1", "Kdr", "Emcn"
  ),
  group.by = "celltype"
) +
  RotatedAxis()
DimPlot(
  seu_m,
  group.by = "celltype",
  label    = TRUE,
  repel    = TRUE
) +
  ggtitle("Panel A – Cell-Type Annotation of Mouse Limb Growth Plate") +
  theme(legend.position = "right")
markers.to.plot.m <- c(
  "Sox9", "Col2a1", "Acan", "Prrx1", "Col1a1", "Twist1",
  "Pecam1", "Kdr", "Emcn", "Fli1",
  "Pitx1", "Wnt7a", "Tfap2a", "Grhl2", "Lama5"
)

dp_m <- DotPlot(
  seu_m,
  features = markers.to.plot.m,
  group.by = "celltype"
)

dp_data_m <- dp_m$data

ggplot(dp_data_m, aes(x = id, y = features.plot)) +
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
    title = "Mouse – Dotplot of Network Genes",
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
markers.to.plot.m2 <- c(
  # 1) Sender: ECM + ligand
  "Lama5", "Wnt7a",
  # 2) Boundary
  "Tfap2a", "Grhl2",
  # 3) Receiver: transcriptional targets
  "Pitx1", "Fli1",
  # 4) Identity – chondrocytes first (receiver), then perichondrium (sender base)
  "Sox9", "Col2a1", "Acan",
  "Prrx1", "Col1a1", "Twist1",
  # 5) Vascular reference
  "Pecam1", "Kdr", "Emcn"
)

dp_m2 <- DotPlot(
  seu_m,
  features = markers.to.plot.m2,
  group.by = "celltype"
)

dp_data_m2 <- dp_m2$data

dp_data_m2$id <- factor(
  dp_data_m2$id,
  levels = c("Perichondrium", "Chondrocytes", "Endothelium")
)

dp_data_m2$features.plot <- factor(
  dp_data_m2$features.plot,
  levels = rev(markers.to.plot.m2)
)

ggplot(dp_data_m2, aes(x = id, y = features.plot)) +
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
    title = "Mouse – Perichondrial Sender / Chondrocyte Receiver / Vascular Reference",
    x     = NULL,
    y     = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title       = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.text.x      = element_text(size = 14, angle = 45, hjust = 1),
    axis.text.y      = element_text(size = 12, face = "italic"),
    legend.title     = element_text(size = 14),
    legend.text      = element_text(size = 12),
    panel.grid.major = element_line(color = "grey90")
  )
# 0. Expression matrix for key niche genes
genes_use_m <- c(
  "Lama5", "Wnt7a",
  "Tfap2a", "Grhl2",
  "Pitx1", "Fli1"
)

expr_data_m <- t(as.matrix(
  GetAssayData(seu_m, assay = "RNA", slot = "data")[genes_use_m, , drop = FALSE]
))

stopifnot(nrow(expr_data_m) == nrow(seu_m@meta.data))

celltype_m <- as.character(seu_m$celltype)

# 1. Hierarchy
seu_m$Niche_Architecture <- vapply(
  X = seq_len(nrow(expr_data_m)),
  FUN = function(i) {
    x <- expr_data_m[i, ]
    
    # --- GROUP A: ARCHITECTS ---
    if (x["Lama5"] > 0 && x["Wnt7a"] > 0) {
      return("Niche Hub (Lama5+/Wnt7a+)")
    }
    if (x["Lama5"] > 0) {
      return("Matrix Architect (Lama5+)")
    }
    if (x["Wnt7a"] > 0) {
      return("Niche Signal (Wnt7a+)")
    }
    if (x["Tfap2a"] > 0 || x["Grhl2"] > 0) {
      return("Boundary (Tfap2a+/Grhl2+)")
    }
    
    # --- GROUP B: TARGETS ---
    if (x["Pitx1"] > 0 && x["Fli1"] > 0) {
      return("Target: Dual Positive (Pitx1+/Fli1+)")
    }
    if (x["Pitx1"] > 0) {
      return("Target: Pitx1+")
    }
    if (x["Fli1"] > 0) {
      return("Target: Fli1+")
    }
    
    # --- BACKGROUND: coarse celltype
    return(celltype_m[i])
  },
  FUN.VALUE = character(1)
)
plot_order_m <- c(
  "Niche Hub (Lama5+/Wnt7a+)",
  "Niche Signal (Wnt7a+)",
  "Matrix Architect (Lama5+)",
  "Boundary (Tfap2a+/Grhl2+)",
  "Target: Dual Positive (Pitx1+/Fli1+)",
  "Target: Pitx1+",
  "Target: Fli1+",
  "Chondrocytes",
  "Perichondrium",
  "Endothelium"
)

niche_colors_m <- c(
  "Niche Hub (Lama5+/Wnt7a+)"           = "green",
  "Niche Signal (Wnt7a+)"               = "#4B0092",
  "Matrix Architect (Lama5+)"           = "#332288",
  "Boundary (Tfap2a+/Grhl2+)"           = "#88CCEE",
  "Target: Dual Positive (Pitx1+/Fli1+)"= "#D55E00",
  "Target: Pitx1+"                      = "#E69F00",
  "Target: Fli1+"                       = "#CC79A7",
  "Chondrocytes"                        = "#F0E442",
  "Perichondrium"                       = "#DDDDDD",
  "Endothelium"                         = "grey70",
  "Other"                               = "grey98"
)

p_m <- DimPlot(
  seu_m,
  group.by = "Niche_Architecture",
  cols     = niche_colors_m,
  order    = plot_order_m,
  pt.size  = 1.2
) +
  ggtitle("Mouse Signaling Niche") +
  theme_void() +
  theme(
    plot.title      = element_text(size = 18, face = "bold", hjust = 0.5),
    legend.position = "right",
    legend.title    = element_blank(),
    legend.text     = element_text(size = 18)
  )
umap_mat_m <- Embeddings(seu_m, "umap")
xr_m <- range(umap_mat_m[, 1])
yr_m <- range(umap_mat_m[, 2])

dx_m <- diff(xr_m) * 0.20
dy_m <- diff(yr_m) * 0.20

p_m_with_axes <- p_m +
  geom_segment(
    aes(x = xr_m[1], y = yr_m[1], xend = xr_m[1] + dx_m, yend = yr_m[1]),
    arrow       = arrow(length = unit(0.2, "cm"), type = "closed"),
    linewidth   = 0.8,
    inherit.aes = FALSE
  ) +
  geom_segment(
    aes(x = xr_m[1], y = yr_m[1], xend = xr_m[1], yend = yr_m[1] + dy_m),
    arrow       = arrow(length = unit(0.2, "cm"), type = "closed"),
    linewidth   = 0.8,
    inherit.aes = FALSE
  ) +
  annotate(
    "text",
    x     = xr_m[1],
    y     = yr_m[1] - 0.02 * diff(yr_m),
    label = "UMAP1",
    hjust = 0,
    vjust = 1,
    size  = 6
  ) +
  annotate(
    "text",
    x     = xr_m[1] - 0.02 * diff(xr_m),
    y     = yr_m[1],
    label = "UMAP2",
    hjust = 0,
    vjust = 0,
    angle = 90,
    size  = 6
  )

p_m_with_axes
celltype_colors_m <- c(
  "Chondrocytes"  = "#F0E442",
  "Perichondrium" = "#DDDDDD",
  "Endothelium"   = "grey70"
)

ct_m <- DimPlot(
  seu_m,
  group.by   = "celltype",
  cols       = celltype_colors_m,
  label      = FALSE,
  label.size = 6
) +
  ggtitle("Mouse Cell Type Annotation") +
  theme_void() +
  theme(
    plot.title      = element_text(size = 18, face = "bold", hjust = 0.5),
    legend.position = "right",
    legend.title    = element_blank(),
    legend.text     = element_text(size = 18),
  )

umap_mat_m <- Embeddings(seu_m, "umap")
xr_m <- range(umap_mat_m[, 1])
yr_m <- range(umap_mat_m[, 2])

dx_m <- diff(xr_m) * 0.20
dy_m <- diff(yr_m) * 0.20

ct_m_with_axes <- ct_m +
  geom_segment(
    aes(x = xr_m[1], y = yr_m[1], xend = xr_m[1] + dx_m, yend = yr_m[1]),
    arrow       = arrow(length = unit(0.25, "cm"), type = "closed"),
    linewidth   = 0.9,
    inherit.aes = FALSE
  ) +
  geom_segment(
    aes(x = xr_m[1], y = yr_m[1], xend = xr_m[1], yend = yr_m[1] + dy_m),
    arrow       = arrow(length = unit(0.25, "cm"), type = "closed"),
    linewidth   = 0.9,
    inherit.aes = FALSE
  ) +
  annotate(
    "text",
    x     = xr_m[1],
    y     = yr_m[1] - 0.02 * diff(yr_m),
    label = "UMAP1",
    hjust = 0,
    vjust = 1,
    size  = 6
  ) +
  annotate(
    "text",
    x     = xr_m[1] - 0.02 * diff(xr_m),
    y     = yr_m[1],
    label = "UMAP2",
    hjust = 0,
    vjust = 0,
    angle = 90,
    size  = 6
  )

ct_m_with_axes
