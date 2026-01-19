library(ggplot2)
library(dplyr)
library(ggrepel)
library(readxl)

# ------------------- General Volcano Plot Function ------------------- #

create_volcano_plot <- function(file_path, condition_label, marker_genes, colors, title_text = NULL) {
  
  # Load and clean table
  deg_table <- read_xlsx(file_path)
  deg_table_clean <- deg_table %>%
    filter(!is.na(log2FoldChange), !is.na(padj))
  deg_table_clean$log2FoldChange <- -deg_table_clean$log2FoldChange
  
  # Thresholds
  pval_threshold <- 0.05
  logfc_threshold <- 1
  
  # Categorize by significance and direction
  deg_table_clean <- deg_table_clean %>%
    mutate(
      Significance = case_when(
        padj < pval_threshold & abs(log2FoldChange) > logfc_threshold ~ "Significant",
        TRUE ~ "Not significant"
      ),
      Direction = case_when(
        padj < pval_threshold & log2FoldChange > logfc_threshold ~ paste0(condition_label, " up"),
        padj < pval_threshold & log2FoldChange < -logfc_threshold ~ paste0(condition_label, " down"),
        TRUE ~ "Not significant"
      ),
      Color = Direction,
      Marker = ifelse(GeneSymbol %in% marker_genes, GeneSymbol, NA)
    )
  
  # Filter extreme outliers
  filtered_deg_table <- deg_table_clean %>%
    filter(between(log2FoldChange, -10, 10))
  
  # Plot
  # Plot
  volcano_plot <- ggplot(filtered_deg_table, aes(x = log2FoldChange, y = -log10(padj), color = Color)) +
    geom_point(alpha = 0.85, size = 3) +
    scale_color_manual(values = colors) +
    geom_vline(xintercept = c(-logfc_threshold, logfc_threshold), 
               linetype = "dotted", color = "grey40", linewidth = 0.6) +
    geom_hline(yintercept = -log10(pval_threshold), 
               linetype = "dotted", color = "grey40", linewidth = 0.6) +
    geom_point(data = filtered_deg_table %>% filter(!is.na(Marker)), 
               size = 4.2, color = "black", shape = 21, fill = "white", stroke = 1.2) +
    geom_label_repel(
      data = filtered_deg_table %>% filter(!is.na(Marker)),
      aes(label = Marker),
      color = "black", 
      size = 6,
      fill = "white",
      fontface = "bold.italic",
      box.padding = 0.3, 
      point.padding = 0.3,
      max.overlaps = 20
    ) +
    labs(
      title = title_text,
      x = expression(log[2]~"Fold Change"),
      y = expression(-log[10]~"(p-adj.)"),
      color = NULL
    ) +
    theme_classic(base_size = 16) +
    theme(
      legend.position = "none",
      plot.title      = element_text(size = 18, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 16, face = "bold"),
      axis.text = element_text(size = 14),
      axis.line = element_line(size = 0.8),
      panel.grid = element_blank()
    )
  
  return(volcano_plot)
}



# ------------------- Plot CDH1 ------------------- #

DEG_path <- "...differentialExpression_ThielRNA.xlsx"
DEG_markers <- c("FLI1", "WNT7A")
DEG_colors <- c("LAMA5 KO up" = "#FF6666", 
                "LAMA5 KO down" = "#66A3FF", 
                "Not significant" = "grey80")

volcano_DEG <- create_volcano_plot(
  file_path = DEG_path,
  condition_label = "LAMA5 KO",
  marker_genes = DEG_markers,
  colors = DEG_colors,
  title_text = "LAMA5 KO vs. WT"
)

print(volcano_DEG)
