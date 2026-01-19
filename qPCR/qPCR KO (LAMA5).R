library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(broom)
library(writexl)

# Load and clean LAMA5 data
file_path <- "...LAMA5 KO validation all celltypes_QuantStudio 12K Flex_export.xlsx"

data <- read_excel(file_path, sheet = 3, col_names = TRUE) %>%
  mutate(
    CT = as.character(CT),
    CT = str_replace(CT, ",", "."),
    CT = ifelse(CT == "Undetermined", "50", CT),
    CT = as.numeric(CT)
  )

# Outlier removal function
remove_outliers <- function(x) {
  Q1 <- quantile(x, 0.25, na.rm = TRUE)
  Q3 <- quantile(x, 0.75, na.rm = TRUE)
  IQR_val <- Q3 - Q1
  lower <- Q1 - 1.5 * IQR_val
  upper <- Q3 + 1.5 * IQR_val
  x[x >= lower & x <= upper]
}

# Updated ΔΔCt processing with outlier removal
process_lama5_data <- function(data, target_gene = "LAMA5", housekeeping_gene = "RPLP0") {
  
  # Add Lineage column
  data <- data %>%
    mutate(Lineage = case_when(
      str_detect(Sample, "^ch") ~ "ch",
      str_detect(Sample, "^o") ~ "o",
      TRUE                     ~ ""
    ))
  
  # Keep only target + housekeeping
  data_filtered <- data %>%
    filter(Target %in% c(target_gene, housekeeping_gene)) %>%
    group_by(Lineage, Sample, Target) %>%
    mutate(rep = row_number()) %>%
    ungroup()
  
  # 🧹 Outlier removal
  data_cleaned <- data_filtered %>%
    group_by(Sample, Target) %>%
    mutate(CT_clean = list(remove_outliers(CT))) %>%
    ungroup() %>%
    rowwise() %>%
    mutate(CT_filtered = ifelse(length(CT_clean) > 0, list(CT_clean), list(NA))) %>%
    unnest(CT_filtered) %>%
    select(Sample, Target, CT = CT_filtered, Lineage)
  
  # 🚩 Flag samples that still contain CT = 50 (ND) AFTER outlier removal
  nd_flag <- data_cleaned %>%
    filter(CT == 50) %>%
    distinct(Sample) %>%
    mutate(Contains_ND = TRUE)
  
  # Pivot to wide format
  data_wide <- data_cleaned %>%
    group_by(Lineage, Sample, Target) %>%
    mutate(rep = row_number()) %>%
    ungroup() %>%
    pivot_wider(
      id_cols = c(Lineage, Sample, rep),
      names_from = Target,
      values_from = CT
    ) %>%
    rename(Gene_Ct = !!target_gene, HK_Ct = !!housekeeping_gene) %>%
    filter(!is.na(Gene_Ct) & !is.na(HK_Ct))
  
  # Mean per biological replicate
  data_biorep <- data_wide %>%
    group_by(Lineage, Sample) %>%
    summarise(
      Gene_Ct_mean = mean(Gene_Ct, na.rm = TRUE),
      HK_Ct_mean   = mean(HK_Ct, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      Delta_Ct = Gene_Ct_mean - HK_Ct_mean,
      Genotype = if_else(str_detect(Sample, "KO"), "KO", "WT"),
      Gene = target_gene
    ) %>%
    # ⬅️ Join ND flag
    left_join(nd_flag, by = "Sample") %>%
    mutate(Contains_ND = ifelse(is.na(Contains_ND), FALSE, Contains_ND))
  
  # Reference ΔCt from WT
  reference_dCt <- data_biorep %>%
    filter(Genotype == "WT") %>%
    group_by(Lineage) %>%
    summarise(mean_Delta_Ct = mean(Delta_Ct, na.rm = TRUE), .groups = "drop")
  
  # ΔΔCt and RQ
  final_data <- data_biorep %>%
    left_join(reference_dCt, by = "Lineage") %>%
    mutate(
      DeltaDelta_Ct = Delta_Ct - mean_Delta_Ct,
      Normalized_RQ = 2^(-DeltaDelta_Ct)
    )
  
  # T-test: KO vs WT per lineage
  ttest_results <- final_data %>%
    group_by(Lineage) %>%
    do(tidy(t.test(DeltaDelta_Ct ~ Genotype, data = .))) %>%
    mutate(Gene = target_gene) %>%
    ungroup()
  
  return(list(
    results = final_data,
    ttest   = ttest_results
  ))
}

# Run analysis
lama5_result <- process_lama5_data(data, "LAMA5", "RPLP0")

# Export
output_path <- "...normalized_RQ_LAMA5_Enhanced.xlsx"

write_xlsx(list(
  normalized_RQ = lama5_result$results,
  ttest_summary = lama5_result$ttest
), path = output_path)

message("Enhanced LAMA5 ΔΔCt analysis with outlier removal saved to: ", output_path)
