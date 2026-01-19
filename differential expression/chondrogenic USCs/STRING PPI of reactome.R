library(ReactomePA)
library(org.Hs.eg.db)
library(STRINGdb)
library(tidyverse)
library(ggraph)
library(igraph)
library(reactome.db)
library(AnnotationDbi)

#─── 1) Get Entrez IDs for Reactome pathway R‑HSA‑446728 ────────────────────
pathway_id <- "R-HSA-446728"
reactome_anno <- AnnotationDbi::select(
  x       = reactome.db,
  keys    = pathway_id,
  keytype = "PATHID",
  columns = "ENTREZID"
)
genes_entrez <- unique(reactome_anno$ENTREZID)

#─── 2) Read DE table and map ENSEMBL → ENTREZID & SYMBOL ───────────────────
deg_table <- readxl::read_excel(
  "...differentialExpression_Thiel_Chond.xlsx"
) %>%
  mutate(
    ENSEMBL = sub("\\..*", "", ENSEMBL),   # Remove version suffix from ENSEMBL IDs
    log2FC  = -log2FoldChange              # Reverse direction to KO vs WT
  )

deg_entrez <- AnnotationDbi::select(
  x       = org.Hs.eg.db,
  keys    = deg_table$ENSEMBL,
  keytype = "ENSEMBL",
  columns = c("ENTREZID", "SYMBOL")
) %>%
  rename(gene_symbol = SYMBOL) %>%         # rename here to avoid confusion
  inner_join(deg_table, by = "ENSEMBL")


# Subset to pathway genes and make named vector of log2FC
deg_pathway <- deg_entrez %>%
  filter(ENTREZID %in% genes_entrez)

log2fc_vec  <- setNames(deg_pathway$log2FC, deg_pathway$ENTREZID)

#─── 3) Fetch high‑confidence PPIs among these genes from STRING ───────────
string_db <- STRINGdb$new(
  version         = "11",
  species         = 9606,
  score_threshold = 700,
  input_directory = ""
)

# map Entrez → STRING IDs
mapped <- string_db$map(
  data.frame(ENTREZID = genes_entrez),
  "ENTREZID",
  removeUnmappedRows = TRUE
)

# get interactions where both partners are in our list
ppi_all <- string_db$get_interactions(mapped$STRING_id)

# 3a) Inspect the columns
print(head(ppi_all))
print(colnames(ppi_all))

# 3b) Suppose we see columns called "from" and "to" (or stringId_A/stringId_B),
#     adjust the filter accordingly. For example, if the names are "from"/"to":

ppi_sub <- ppi_all %>%
  filter(from %in% mapped$STRING_id,
         to   %in% mapped$STRING_id)



#─── 4) Build igraph object & annotate with log2FC + gene symbols ──────────
# Build igraph directly from the filtered edges:
# Directly get the subgraph of your IDs at score ≥ 700
library(igraph)

# Build igraph from STRING PPI data frame
g <- igraph::graph_from_data_frame(
  d = ppi_sub,
  directed = FALSE
)


# Merge mapped (STRING → Entrez) with deg_pathway to get annotation per node
node_annot <- mapped %>%
  inner_join(deg_pathway, by = "ENTREZID") %>%
  select(STRING_id, gene_symbol, log2FC)


# Assign attributes to nodes
V(g)$log2FC <- node_annot$log2FC[match(V(g)$name, node_annot$STRING_id)]
V(g)$label  <- node_annot$gene_symbol[match(V(g)$name, node_annot$STRING_id)]

library(ggraph)
library(igraph)
library(scales)
library(dplyr)

# ── 1) Remove nodes with missing or small log2FC ─────────────────────
g_clean <- igraph::delete_vertices(g, V(g)[is.na(log2FC) | abs(log2FC) < 1])

# ── 2) Keep only the largest connected component ─────────────────────
components <- igraph::components(g_clean)
main_component_nodes <- V(g_clean)$name[components$membership == which.max(components$csize)]
g_main <- igraph::induced_subgraph(g_clean, vids = main_component_nodes)

# ── 3) Annotate edge weights (STRING confidence scores) ──────────────

# Extract edge list from the igraph object
edges_main_clean <- as_data_frame(g_main, what = "edges") %>%
  mutate(
    from_u = pmin(from, to),  # undirected edge key: alphabetically first
    to_u   = pmax(from, to),
    key    = paste(from_u, to_u, sep = "_")
  )

# Do the same normalization of keys for the PPI data
ppi_sub_clean <- ppi_sub %>%
  mutate(
    from_u = pmin(from, to),
    to_u   = pmax(from, to),
    key    = paste(from_u, to_u, sep = "_")
  )

# Create named vector of combined_scores
score_vec <- setNames(ppi_sub_clean$combined_score, ppi_sub_clean$key)

# Map scores back to edges_main_clean
edges_main_clean$combined_score <- score_vec[edges_main_clean$key]

# Check for NAs
sum(is.na(edges_main_clean$combined_score))  # Ideally 0

# Assign weights
E(g_main)$weight <- edges_main_clean$combined_score

# Scale weights
E(g_main)$weight_scaled <- scales::rescale(E(g_main)$weight, to = c(0.2, 2))


# Raw values
node_sizes_raw <- abs(V(g_main)$log2FC)

# Rescale to a reasonable size range, e.g., 3 to 8 (this controls dot sizes)
node_sizes_scaled <- scales::rescale(node_sizes_raw, to = c(4, 10))

# Assign back to graph nodes
V(g_main)$size_scaled <- node_sizes_scaled

# ── 4) Plot the graph ────────────────────────────────────────────────
library(ggraph)
library(ggrepel)
# Split layout data
layout_df <- create_layout(g_main, layout = "fr")

highlight_gene <- "CDH1"

label_cdh1 <- layout_df %>% filter(label == highlight_gene)
label_others <- layout_df %>% filter(label != highlight_gene)

ggraph(layout_df) +
  geom_edge_link(aes(width = weight_scaled), alpha = 0.3, color = "gray80") +
  geom_node_point(aes(color = log2FC, size = size_scaled)) +
  
  # Label for CDH1
  ggrepel::geom_label_repel(
    data = label_cdh1,
    aes(x = x, y = y, label = label),
    fill        = "white",
    color       = "black",
    size        = 6,  # bigger label
    fontface    = "bold.italic",
    box.padding = 0.3,
    point.padding = 0.3,
    max.overlaps = Inf
  ) +
  
  # Labels for others (smaller, lighter)
  ggrepel::geom_label_repel(
    data = label_others,
    aes(x = x, y = y, label = label),
    fill        = "white",
    color       = "black",
    size        = 4,  # smaller
    fontface    = "italic",
    box.padding = 0.2,
    point.padding = 0.2,
    max.overlaps = 10
  ) +
  
  scale_edge_width(range = c(0.2, 2), guide = "none") +
  scale_size_continuous(range = c(4, 10), guide = "none") +
  scale_color_gradient2(
    low = "navy",
    mid = "white",
    high = "firebrick3",
    midpoint = 0,
    name = "log2FC"
  ) +
  theme_void() +
  labs(
    title    = "Cell junction organization PPI"
  ) +
  theme(
    plot.title    = element_text(size = 18, face = "bold", hjust = 0.5),
    legend.title  = element_text(size = 14, face = "bold"),
    legend.text   = element_text(size = 12)
  )

