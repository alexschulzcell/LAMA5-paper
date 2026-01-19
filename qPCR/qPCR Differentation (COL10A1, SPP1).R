library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(broom)
library(writexl)

# 🧹 Outlier removal
remove_outliers <- function(x) {
  Q1 <- quantile(x, 0.25, na.rm = TRUE)
  Q3 <- quantile(x, 0.75, na.rm = TRUE)
  IQR <- Q3 - Q1
  lower <- Q1 - 1.5 * IQR
  upper <- Q3 + 1.5 * IQR
  x[x >= lower & x <= upper]
}

# 📦 Load and clean qPCR data
load_qpcr <- function(path, sheet) {
  read_excel(path, sheet = sheet) %>%
    mutate(
      CT = as.character(CT),
      CT = str_replace_all(CT, ",", "."),
      CT = ifelse(CT == "Undetermined", "50", CT),
      CT = as.numeric(CT)
    )
}

# 🔬 Enhanced ΔΔCt processing for differentiation comparison
process_diff_comparison <- function(df, gene, hk = "RPLP0") {
  df2 <- df %>%
    mutate(
      Lineage = case_when(
        str_detect(Sample, "^ch") ~ "ch",
        str_detect(Sample, "^o")  ~ "o",
        TRUE                      ~ ""   # undifferentiated
      ),
      Genotype = case_when(
        str_detect(Sample, "WT") ~ "WT",
        str_detect(Sample, "KO") ~ "KO",
        TRUE                     ~ NA_character_
      )
    ) %>%
    filter(Target %in% c(gene, hk))
  
  # Outlier removal
  df2_clean <- df2 %>%
    group_by(Sample, Target) %>%
    mutate(CT_clean = list(remove_outliers(CT))) %>%
    ungroup() %>%
    rowwise() %>%
    mutate(CT_filtered = ifelse(length(CT_clean) > 0, list(CT_clean), list(NA))) %>%
    unnest(CT_filtered) %>%
    select(Sample, Target, CT = CT_filtered) %>%
    left_join(df2 %>% select(Sample, Target, Lineage, Genotype) %>% distinct(), 
              by = c("Sample", "Target"))
  
  # Re-flag samples that still contain ND (CT = 50)
  nd_flags <- df2_clean %>%
    filter(CT == 50) %>%
    distinct(Sample) %>%
    mutate(Contains_ND = TRUE)
  
  # Continue as before
  df2_clean <- df2_clean %>%
    group_by(Lineage, Genotype, Sample, Target) %>%
    mutate(rep = row_number()) %>%
    ungroup() %>%
    pivot_wider(
      id_cols     = c(Lineage, Genotype, Sample, rep),
      names_from  = Target,
      values_from = CT
    ) %>%
    rename(
      Gene_Ct = all_of(gene),
      HK_Ct   = all_of(hk)
    ) %>%
    filter(!is.na(Gene_Ct), !is.na(HK_Ct))
  
  # Calculate means
  df_summary <- df2_clean %>%
    group_by(Lineage, Genotype, Sample) %>%
    summarise(
      Gene_Ct_mean = mean(Gene_Ct, na.rm = TRUE),
      HK_Ct_mean   = mean(HK_Ct, na.rm = TRUE),
      .groups      = "drop"
    ) %>%
    mutate(Delta_Ct = Gene_Ct_mean - HK_Ct_mean) %>%
    left_join(nd_flags, by = "Sample") %>%
    mutate(Contains_ND = ifelse(is.na(Contains_ND), FALSE, Contains_ND))
  
  # Calibrator: based on lineage
  cal_lineage <- if (gene == "SPP1") "o" else "ch"
  calibrator_delta <- df_summary %>%
    filter(Lineage == cal_lineage) %>%
    pull(Delta_Ct) %>%
    mean(na.rm = TRUE)
  
  # Final result
  df_final <- df_summary %>%
    mutate(
      DeltaDelta_Ct = Delta_Ct - calibrator_delta,
      Normalized_RQ = 2^(-DeltaDelta_Ct),
      Gene          = gene
    ) %>%
    select(Gene, Genotype, Sample, Lineage, Gene_Ct_mean, HK_Ct_mean,
           Delta_Ct, DeltaDelta_Ct, Normalized_RQ, Contains_ND)
  
  return(df_final)
}


# 🧪 T-tests: Undiff vs Diff per genotype
diff_ttest <- function(df_all) {
  df_all %>%
    filter(!is.na(Genotype)) %>%
    group_by(Gene, Genotype) %>%
    do({
      sub <- .
      if (n_distinct(sub$Lineage) >= 2) {
        tidy(t.test(DeltaDelta_Ct ~ Lineage, data = sub))
      } else {
        tibble(
          statistic   = NA_real_,
          p.value     = NA_real_,
          conf.low    = NA_real_,
          conf.high   = NA_real_,
          method      = "t-test",
          alternative = "two.sided"
        )
      }
    }) %>%
    ungroup()
}


# Define input files
gene_files <- list(
  SPP1 = list(
    path  = "...RNA Seq Validation osteos_QuantStudio 12K Flex_export.xlsx",
    sheet = 2
  ),
  COL10A1 = list(
    path  = "...Alex qPCR RNA Seq Validierung_QuantStudio 12K Flex_export.xlsx",
    sheet = 5
  )
)

# Process each gene
gene_results <- lapply(names(gene_files), function(gene) {
  info <- gene_files[[gene]]
  df_raw <- load_qpcr(info$path, info$sheet)
  process_diff_comparison(df_raw, gene)
})
names(gene_results) <- names(gene_files)

# Combine all for stats
all_results <- bind_rows(gene_results)
ttest_results <- diff_ttest(all_results)

# Export
output_path <- "...RQ_diff_vs_undiff_SPP1_COL10A1_Enhanced.xlsx"

write_xlsx(
  c(
    setNames(gene_results, paste0(names(gene_results), "_RQ")),
    list(TTest_Results = ttest_results)
  ),
  path = output_path
)

message("✅ Enhanced SPP1/COL10A1 ΔΔCt diff analysis exported to: ", output_path)


