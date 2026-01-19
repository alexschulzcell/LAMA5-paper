library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(broom)
library(writexl)

# Outlier removal function (from first script)
remove_outliers <- function(x) {
  Q1 <- quantile(x, 0.25, na.rm = TRUE)
  Q3 <- quantile(x, 0.75, na.rm = TRUE)
  IQR_val <- Q3 - Q1
  lower_bound <- Q1 - 1.5 * IQR_val
  upper_bound <- Q3 + 1.5 * IQR_val
  x[x >= lower_bound & x <= upper_bound]
}

# Unified function
process_gene_combined <- function(data, gene_name, housekeeping_gene = "RPLP0", cell_type = "Undifferentiated") {
  
  # Rename and filter by cell type
  data <- data %>%
    rename(Sample = ...1, Target = ...2, CT = ...3) %>%
    filter(
      (cell_type == "Undifferentiated" & !str_starts(Sample, "ch") & !str_starts(Sample, "o")) |
        (cell_type == "Chondrogenic" & str_starts(Sample, "ch")) |
        (cell_type == "Osteogenic" & str_starts(Sample, "o"))
    )
  
  # Clean CT values
  data <- data %>%
    mutate(
      CT = as.character(CT),
      CT = str_replace_all(CT, ",", "."),
      CT = ifelse(CT == "Undetermined", "50", CT),
      CT = as.numeric(CT)
    )
  
  # Outlier removal per Sample and Target
  data_cleaned <- data %>%
    group_by(Sample, Target) %>%
    mutate(
      cleaned_CT = list(remove_outliers(CT))
    ) %>%
    ungroup() %>%
    rowwise() %>%
    mutate(CT_filtered = ifelse(length(cleaned_CT) > 0, list(cleaned_CT), list(NA))) %>%
    unnest(CT_filtered) %>%
    select(Sample, Target, CT = CT_filtered)
  
  # Add Lineage column
  data_cleaned <- data_cleaned %>%
    mutate(Lineage = case_when(
      cell_type == "Osteogenic" ~ "o",
      cell_type == "Chondrogenic" ~ "ch",
      TRUE ~ ""
    ))
  # Re-flag samples with ND (CT = 50) that survived filtering
  nd_flags <- data_cleaned %>%
    filter(CT == 50) %>%
    distinct(Sample) %>%
    mutate(Contains_ND = TRUE)
  
  # Keep only gene of interest + housekeeping
  data_filtered <- data_cleaned %>%
    filter(Target %in% c(gene_name, housekeeping_gene)) %>%
    group_by(Lineage, Sample, Target) %>%
    mutate(rep = row_number()) %>%
    ungroup()
  
  # Pivot to wide
  data_wide <- data_filtered %>%
    pivot_wider(
      id_cols = c(Lineage, Sample, rep),
      names_from = Target,
      values_from = CT
    ) %>%
    rename(Gene_Ct = !!gene_name, HK_Ct = !!housekeeping_gene) %>%
    filter(!is.na(Gene_Ct) & !is.na(HK_Ct))
  
  # Calculate mean CTs
  data_biorep <- data_wide %>%
    group_by(Lineage, Sample) %>%
    summarise(
      Gene_Ct_mean = mean(Gene_Ct, na.rm = TRUE),
      HK_Ct_mean = mean(HK_Ct, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      Delta_Ct = Gene_Ct_mean - HK_Ct_mean,
      Genotype = if_else(str_detect(Sample, "KO"), "KO", "WT"),
      Cell_Type = cell_type,
      Gene = gene_name
    )
  
  # Get WT reference for ΔΔCt
  reference_dCt <- data_biorep %>%
    filter(Genotype == "WT") %>%
    group_by(Lineage) %>%
    summarise(mean_Delta_Ct = mean(Delta_Ct, na.rm = TRUE), .groups = "drop")
  
  # Calculate ΔΔCt and RQ
  data_final <- data_biorep %>%
    left_join(reference_dCt, by = "Lineage") %>%
    mutate(
      DeltaDelta_Ct = Delta_Ct - mean_Delta_Ct,
      RQ = 2^(-DeltaDelta_Ct)
    )
  
  # Normalize RQ per gene logic (like first script)
  data_final <- data_final %>%
    group_by(Lineage) %>%
    mutate(
      reference_group = case_when(
        Gene == "FLI1" ~ "KO",
        Gene == "WNT7A" ~ "WT",
        TRUE ~ NA_character_
      ),
      norm_factor = if (!is.na(reference_group[1])) {
        mean(RQ[Genotype == reference_group[1]], na.rm = TRUE)
      } else {
        NA_real_
      },
      Normalized_RQ = ifelse(!is.na(norm_factor), RQ / norm_factor, NA_real_)
    ) %>%
    ungroup()
  
  # T-test: KO vs WT
  ttest_results <- data_final %>%
    group_by(Lineage) %>%
    do(tidy(t.test(DeltaDelta_Ct ~ Genotype, data = .))) %>%
    mutate(Gene = gene_name, Cell_Type = cell_type) %>%
    ungroup()
  # Add ND flag to final results
  data_final <- data_final %>%
    left_join(nd_flags, by = "Sample") %>%
    mutate(Contains_ND = ifelse(is.na(Contains_ND), FALSE, Contains_ND))
  
  
  return(list(
    results = data_final,
    ttest = ttest_results
  ))
}



# File path
file_path <- "...Alex qPCR RNA Seq Validierung_QuantStudio 12K Flex_export.xlsx"

# Load sheets
data_list <- list(
  FLI1_Undiff = read_excel(file_path, sheet = 6, col_names = FALSE),
  FLI1_Chondro = read_excel(file_path, sheet = 6, col_names = FALSE),
  WNT7A_Undiff = read_excel(file_path, sheet = 7, col_names = FALSE),
  WNT7A_Chondro = read_excel(file_path, sheet = 7, col_names = FALSE),
  FLI1_Osteo = read_excel(file_path, sheet = 10, col_names = FALSE),
  WNT7A_Osteo = read_excel(file_path, sheet = 11, col_names = FALSE)
)

# Run processing
results_list <- list()
ttest_list <- list()

for (name in names(data_list)) {
  gene <- ifelse(str_detect(name, "FLI1"), "FLI1", "WNT7A")
  cell_type <- case_when(
    str_detect(name, "Undiff") ~ "Undifferentiated",
    str_detect(name, "Osteo") ~ "Osteogenic",
    str_detect(name, "Chondro") ~ "Chondrogenic"
  )
  
  res <- process_gene_combined(data_list[[name]], gene, cell_type = cell_type)
  
  sheet_name <- paste(gene, cell_type, sep = "_")
  results_list[[sheet_name]] <- res$results
  ttest_list[[sheet_name]] <- res$ttest
}

# Combine t-tests
ttests_combined <- bind_rows(ttest_list)
results_list[["TTest_Results"]] <- ttests_combined

# Export all to Excel
write_xlsx(results_list, path = "Gene_Expression_Combined_Analysis.xlsx")
