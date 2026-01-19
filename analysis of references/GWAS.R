library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel) # Highly recommended for labels

# 1. LOAD & PREP (Standard)
gwas_file <- "...GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS_EUR.gz"
gwas <- fread(gwas_file)

gwas <- gwas %>%
  mutate(CHR = as.numeric(CHR), POS = as.numeric(POS), P = as.numeric(P), logP = -log10(P)) %>%
  filter(!is.na(CHR), !is.na(POS), !is.na(P), P > 0)

# 2. CALCULATE CUMULATIVE BP ON THE FULL DATASET (Crucial Step)
# This ensures both your background and highlights share the EXACT same coordinate system
chr_lengths <- gwas %>%
  group_by(CHR) %>%
  summarise(chr_max_bp = max(POS)) %>%
  arrange(CHR) %>%
  mutate(chr_offset = cumsum(as.numeric(lag(chr_max_bp, default = 0))))

gwas <- gwas %>%
  left_join(chr_lengths, by = "CHR") %>%
  mutate(BP_cum = POS + chr_offset)

# Define axis ticks
axis_df <- chr_lengths %>%
  mutate(center = chr_offset + chr_max_bp / 2)

# 3. DEFINE HIGHLIGHTS (+/- 250kb)
genes <- data.frame(
  gene = c("WNT7A", "PITX1", "FLI1", "GRHL2", "TFAP2A", "LAMA5"),
  chr  = c(3, 5, 11, 8, 6, 20),
  start = c(13857755, 134363424, 128556430, 102504667, 10393419, 60883011),
  end   = c(13921568, 134370503, 128683162, 102681954, 10419892, 60942368)
)
window <- 250000

gwas$highlight <- "Other"
for (i in seq_len(nrow(genes))) {
  idx <- which(gwas$CHR == genes$chr[i] & 
                 gwas$POS >= (genes$start[i] - window) & 
                 gwas$POS <= (genes$end[i] + window))
  gwas$highlight[idx] <- genes$gene[i]
}

# 4. PRUNING THE BACKGROUND (For speed/aesthetics)
# We keep ALL highlighted SNPs, but only 1 lead SNP per 1MB for the background
gwas_highlights <- gwas %>% filter(highlight != "Other")

gwas_background <- gwas %>%
  filter(highlight == "Other") %>%
  mutate(mb_window = floor(POS / 1e6)) %>%
  group_by(CHR, mb_window) %>%
  slice_max(logP, n = 1, with_ties = FALSE) %>%
  ungroup()

# 5. COLORS & LABELS
gene_colors <- c(
  "WNT7A"  = "#D55E00", "PITX1"  = "#0072B2", "FLI1"   = "#009E73",
  "GRHL2"  = "#CC79A7", "TFAP2A" = "#F0E442", "LAMA5"  = "#E69F00",
  "Other"  = "grey85" # Lighter grey makes colored dots pop
)

# Get top hit per gene for labeling
top_hits <- gwas_highlights %>%
  group_by(highlight) %>%
  slice_max(logP, n = 1, with_ties = FALSE)

# 6. PLOT
ggplot() +
  # Genome-wide Line (Behind points)
  geom_hline(yintercept = -log10(5e-8), linetype = "dashed", color = "grey50", size = 0.3) +
  
  # Background (Pruned)
  geom_point(data = gwas_background, aes(x = BP_cum, y = pmin(logP, 60)), 
             color = "grey85", size = 0.5, alpha = 0.6) +
  
  # Highlights (Full Res)
  geom_point(data = gwas_highlights, aes(x = BP_cum, y = pmin(logP, 60), color = highlight), 
             size = 2.5, alpha = 0.9) +
  
  # Labels
  geom_text_repel(data = top_hits, aes(x = BP_cum, y = pmin(logP, 60), label = highlight, color = highlight),
                  size = 4, fontface = "bold", box.padding = 0.5, nudge_y = 2) +
  
  scale_color_manual(values = gene_colors) +
  scale_x_continuous(breaks = axis_df$center, labels = axis_df$CHR, expand = c(0.01, 0.01)) +
  scale_y_continuous(limits = c(0, 65), expand = c(0, 0)) + # Give room for high peaks
  
  labs(x = "Chromosome", y = expression(-log[10](italic(P))),
       title = "Human Height GWAS (Yengo et al., 2022)",
       subtitle = "The LAMA5-WNT7A-PITX1 network is highly enriched for height-associated variants") +
  
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 8, angle = 0),
    panel.grid.major.y = element_line(color = "grey95"),
    plot.title = element_text(face = "bold")
  )


# 5. COLORS (harmonized with niche UMAP)
gene_colors <- c(
  "WNT7A"  = "#4B0092",  # matches Niche Signal (WNT7A+)
  "PITX1"  = "#E69F00",  # matches Target: PITX1+
  "FLI1"   = "#CC79A7",  # matches Target: FLI1+
  "GRHL2"  = "#88CCEE",  # matches Boundary (TFAP2A+/GRHL2+)
  "TFAP2A" = "#88CCEE",  # same boundary color
  "LAMA5"  = "#332288",  # matches Matrix Architect (LAMA5+)
  "Other"  = "grey90"    # very light background
)

# Top hit per gene for labels
top_hits <- gwas_highlights %>%
  group_by(highlight) %>%
  slice_max(logP, n = 1, with_ties = FALSE) %>%
  ungroup()

# 6. MINIMAL, MATCHY PLOT
p_gwas <- ggplot() +
  # Genome-wide significance line
  geom_hline(
    yintercept = -log10(5e-8),
    linetype   = "dashed",
    color      = "grey70",
    linewidth  = 0.4
  ) +
  
  # Background SNPs (pruned)
  geom_point(
    data  = gwas_background,
    aes(x = BP_cum, y = pmin(logP, 60)),
    color = "grey90",
    size  = 0.4,
    alpha = 0.8
  ) +
  
  # Highlighted loci (full resolution)
  geom_point(
    data = gwas_highlights,
    aes(x = BP_cum, y = pmin(logP, 60), color = highlight),
    size  = 2.3,
    alpha = 0.95
  ) +
  
  # Gene labels (clean, bold)
  geom_text_repel(
    data        = top_hits,
    aes(x = BP_cum, y = pmin(logP, 60), label = highlight, color = highlight),
    size        = 6,
    fontface    = "bold.italic",
    box.padding = 0.5,
    point.padding = 0.2,
    min.segment.length = 0.2,
    seed        = 123
  ) +
  
  scale_color_manual(values = gene_colors) +
  
  scale_x_continuous(
    breaks = axis_df$center,
    labels = axis_df$CHR,
    expand = c(0.01, 0.01)
  ) +
  scale_y_continuous(
    limits = c(0, 65),
    expand = c(0, 0)
  ) +
  
  labs(
    x = "Chromosome",
    y = expression(-log[10](italic(P))),
    title    = "Human Height GWAS (Yengo et al., 2022)",
    subtitle = "The LAMA5–WNT7A–PITX1–FLI1 niche is highly enriched for height-associated variants"
  ) +
  
  theme_classic(base_size = 16) +
  theme(
    # No legend (colors are explained in text/figure panels)
    legend.position = "none",
    
    # X axis text a bit smaller to avoid crowding
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    
    # Light horizontal guide lines, no vertical grid
    panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3),
    panel.grid.minor.y = element_blank(),
    
    plot.title    = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 14, hjust = 0.5)
  )

p_gwas

