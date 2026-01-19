# qPCR pipeline — EXCLUDE one most-outlying BIOREP per Condition (FLI1-based), ANOVA + pairwise Student t-tests vs KO
library(readxl)
library(dplyr)
library(tidyr)
library(writexl)
library(stringr)
library(purrr)
library(broom)

# --------------- USER PARAMETERS ---------------
file_path <- "...LiClrescueFinal2_QuantStudio 12K Flex_export.xlsx"
sheet_to_read <- 2
out_path  <- "...LiClrescueFinal2_RQ_ANOVA_pairwise_vs_KO_withBiorepExclusion.xlsx"

hk_gene <- "RPLP0"            # housekeeping gene
detection_cutoff <- 50        # CT >= this treated as non-detect for statistics
reference_condition <- "KO"   # untreated KO is ΔΔCt reference and pairwise comparator
exclude_one_biorep_per_condition <- TRUE  # toggle the requested exclusion
# ------------------------------------------------

# helper: parse Sample into Condition and Biorep (expects "Condition <space> number" or "CondLabel <num>")
parse_sample_info <- function(sample_vec) {
  tibble(Sample = as.character(sample_vec)) %>%
    mutate(
      Biorep = str_extract(Sample, "\\d+$"),
      Condition = str_trim(str_remove(Sample, "\\s*\\d+$"))
    )
}

# --------------- READ & PREP ---------------
raw <- read_excel(file_path, sheet = sheet_to_read, col_names = TRUE) %>%
  mutate(
    CT_raw = as.character(CT),
    undet = toupper(CT_raw) == "UNDETERMINED",
    CT_str = str_replace(CT_raw, ",", "."),
    CT_num = as.numeric(if_else(undet, NA_character_, CT_str)),
    CT_for_stats = if_else(is.na(CT_num) | CT_num >= detection_cutoff, detection_cutoff, CT_num)
  )

# check required columns
need_cols <- c("Sample", "Target", "CT")
missing_cols <- setdiff(need_cols, names(raw))
if (length(missing_cols) > 0) stop("Missing columns in input: ", paste(missing_cols, collapse = ", "))

# --------------- COLLAPSE TECHNICAL REPLICATES ---------------
collapsed <- raw %>%
  group_by(Sample, Target) %>%
  summarise(
    tech_n = n(),
    tech_nonNA = sum(!is.na(CT_num)),
    tech_any_undet = any(undet, na.rm = TRUE),
    mean_CT = ifelse(tech_nonNA > 0, mean(CT_for_stats[!is.na(CT_num)], na.rm = TRUE), detection_cutoff),
    sd_CT = ifelse(tech_nonNA > 1, sd(CT_for_stats[!is.na(CT_num)], na.rm = TRUE), NA_real_),
    .groups = "drop"
  ) %>%
  mutate(all_ND = tech_nonNA == 0)

# parse samples
sample_info <- parse_sample_info(unique(collapsed$Sample))
collapsed <- collapsed %>% left_join(sample_info, by = "Sample")

# --------------- IDENTIFY one most-outlying BIOREP PER CONDITION (based on FLI1 ΔCt) ---------------
excluded_for_outlier <- character(0)
if (exclude_one_biorep_per_condition) {
  # need mean CTs for FLI1 and HK
  hk_means <- collapsed %>%
    filter(Target == hk_gene) %>%
    select(Sample, mean_CT_hk = mean_CT)
  fli1_means <- collapsed %>%
    filter(Target == "FLI1") %>%
    select(Sample, mean_CT_gene = mean_CT)
  fli1_join <- full_join(fli1_means, hk_means, by = "Sample") %>%
    left_join(sample_info, by = "Sample") %>%
    mutate(Delta_Ct = mean_CT_gene - mean_CT_hk)
  fli1_valid <- fli1_join %>% filter(!is.na(Delta_Ct) & !is.na(Condition))
  if (nrow(fli1_valid) > 0) {
    outlier_list <- fli1_valid %>%
      group_by(Condition) %>%
      mutate(cond_median = median(Delta_Ct, na.rm = TRUE),
             abs_dev = abs(Delta_Ct - cond_median)) %>%
      arrange(desc(abs_dev)) %>%
      slice_head(n = 1) %>%
      ungroup() %>%
      select(Sample, Condition, Biorep, mean_CT_gene, mean_CT_hk, Delta_Ct, abs_dev)
    if (nrow(outlier_list) > 0) {
      message("Identified most outlying biological replicate per Condition (based on FLI1 ΔCt):")
      print(outlier_list)
      excluded_for_outlier <- outlier_list$Sample
    }
  } else {
    message("No valid FLI1 Delta_Ct values found — skipping biorep exclusion.")
  }
}

# remove identified outliers from dataset (if any)
if (length(excluded_for_outlier) > 0) {
  collapsed_filtered <- collapsed %>% filter(!(Sample %in% excluded_for_outlier))
} else {
  collapsed_filtered <- collapsed
}

# --------------- compute per-gene ΔCt and ΔΔCt (reference = untreated KO) ---------------
compute_per_gene <- function(df_collapsed, target_gene, hk_gene = hk_gene, reference_condition = reference_condition) {
  target_df <- df_collapsed %>%
    filter(Target == target_gene) %>%
    select(Sample, Condition, mean_CT_gene = mean_CT, tech_n_gene = tech_n, tech_nonNA_gene = tech_nonNA,
           sd_CT_gene = sd_CT, gene_any_undet = tech_any_undet, gene_all_ND = all_ND)
  hk_df <- df_collapsed %>%
    filter(Target == hk_gene) %>%
    select(Sample, mean_CT_hk = mean_CT, tech_n_hk = tech_n, tech_nonNA_hk = tech_nonNA,
           sd_CT_hk = sd_CT, hk_any_undet = tech_any_undet, hk_all_ND = all_ND)
  joined <- full_join(target_df, hk_df, by = "Sample") %>%
    mutate(
      Gene = target_gene,
      Gene_ND = coalesce(gene_all_ND, FALSE) | coalesce(gene_any_undet, FALSE) | coalesce(mean_CT_gene >= detection_cutoff, FALSE),
      HK_ND   = coalesce(hk_all_ND, FALSE)   | coalesce(hk_any_undet, FALSE)   | coalesce(mean_CT_hk >= detection_cutoff, FALSE),
      Any_ND  = Gene_ND | HK_ND,
      Delta_Ct = mean_CT_gene - mean_CT_hk
    ) %>%
    select(Gene, Sample, Condition,
           mean_CT_gene, mean_CT_hk,
           tech_n_gene, tech_nonNA_gene, sd_CT_gene,
           tech_n_hk, tech_nonNA_hk, sd_CT_hk,
           Gene_ND, HK_ND, Any_ND, Delta_Ct)
  if (nrow(joined) == 0) return(NULL)
  # reference mean Delta_Ct (untreated KO)
  ref_delta_mean <- joined %>%
    filter(Condition == reference_condition) %>%
    summarise(m = mean(Delta_Ct, na.rm = TRUE)) %>%
    pull(m)
  joined <- joined %>%
    mutate(
      DeltaDelta_Ct = ifelse(is.na(Delta_Ct) | is.na(ref_delta_mean), NA_real_, Delta_Ct - ref_delta_mean),
      Normalized_RQ_raw = ifelse(is.na(DeltaDelta_Ct), NA_real_, 2^(-DeltaDelta_Ct))
    )
  # scale so mean of KO raw RQ becomes 1 (if possible)
  ref_mean_rq <- joined %>%
    filter(Condition == reference_condition) %>%
    summarise(m = mean(Normalized_RQ_raw, na.rm = TRUE)) %>%
    pull(m)
  joined <- joined %>%
    mutate(
      Normalized_RQ = ifelse(is.na(Normalized_RQ_raw) | is.na(ref_mean_rq) | ref_mean_rq == 0, NA_real_, Normalized_RQ_raw / ref_mean_rq),
      Mean_Normalized_RQ = Normalized_RQ,
      RQ_ND_flag = Any_ND,
      Mean_Normalized_RQ_for_Prism = if_else(Any_ND, NA_real_, Mean_Normalized_RQ)
    ) %>%
    select(Gene, Sample, Condition,
           mean_CT_gene, mean_CT_hk,
           tech_n_gene, tech_nonNA_gene, sd_CT_gene,
           tech_n_hk, tech_nonNA_hk, sd_CT_hk,
           Gene_ND, HK_ND, Any_ND, RQ_ND_flag,
           Delta_Ct, DeltaDelta_Ct, Normalized_RQ, Mean_Normalized_RQ, Mean_Normalized_RQ_for_Prism)
  return(joined)
}

genes <- sort(unique(collapsed_filtered$Target))
genes <- setdiff(genes, hk_gene)

results_list <- map(genes, ~ compute_per_gene(collapsed_filtered, .x, hk_gene = hk_gene, reference_condition = reference_condition))
names(results_list) <- genes
results_list <- compact(results_list)

combined <- bind_rows(results_list)

# quick check: mean KO RQ per gene should be ~1
ko_check <- combined %>%
  group_by(Gene) %>%
  summarise(mean_RQ_KO = mean(Normalized_RQ[Condition == reference_condition], na.rm = TRUE),
            n_KO = sum(!is.na(Normalized_RQ[Condition == reference_condition]))) %>%
  ungroup()
message("KO mean RQ per gene (should be ~1):"); print(ko_check)

# --------------- SUMMARY stats per gene/condition ---------------
summary_df <- combined %>%
  mutate(ND_flag = Any_ND) %>%
  group_by(Gene, Condition) %>%
  summarise(
    mean_RQ = mean(Normalized_RQ, na.rm = TRUE),
    sd_RQ = sd(Normalized_RQ, na.rm = TRUE),
    n_rep = sum(!is.na(Normalized_RQ)),
    pct_samples_with_any_ND = mean(ND_flag, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Gene, Condition)

# --------------- ANOVA per gene (ΔΔCt ~ Condition) + pairwise Student t-tests vs KO ---------------
anova_list <- list()
pairwise_list <- list()

for (g in sort(unique(combined$Gene))) {
  gene_df <- combined %>% filter(Gene == g)
  usable <- gene_df %>% filter(!Any_ND, !is.na(DeltaDelta_Ct))
  # ANOVA: require at least 3 usable samples and at least two conditions present
  cond_counts <- usable %>% group_by(Condition) %>% summarise(n = n(), .groups = "drop")
  total_usable <- nrow(usable)
  n_conditions_present <- nrow(cond_counts)
  if (total_usable >= 3 && n_conditions_present >= 2) {
    # run aov on DeltaDelta_Ct ~ Condition
    a <- aov(DeltaDelta_Ct ~ Condition, data = usable)
    a_tidy <- broom::tidy(a)
    # extract the Condition row (omnibus)
    cond_row <- a_tidy %>% filter(term == "Condition") %>%
      transmute(Gene = g,
                term = term,
                df1 = df,
                df2 = NA_real_,
                statistic = statistic,
                p.value = p.value)
    if (nrow(cond_row) == 0) {
      s <- summary(a)[[1]]
      if ("Condition" %in% rownames(s)) {
        r <- s["Condition", ]
        cond_row <- tibble(Gene = g, term = "Condition", df1 = as.numeric(r["Df"]), df2 = NA_real_,
                           statistic = as.numeric(r["F value"]), p.value = as.numeric(r["Pr(>F)"]))
      } else {
        cond_row <- tibble(Gene = g, term = "Condition", df1 = NA_real_, df2 = NA_real_, statistic = NA_real_, p.value = NA_real_)
      }
    }
  } else {
    cond_row <- tibble(Gene = g, term = "Condition", df1 = NA_real_, df2 = NA_real_, statistic = NA_real_, p.value = NA_real_)
  }
  anova_list[[length(anova_list) + 1]] <- cond_row
  
  # Pairwise Student t-tests (equal variance) vs KO
  ko_vals <- gene_df %>% filter(Condition == reference_condition, !Any_ND, !is.na(DeltaDelta_Ct)) %>% pull(DeltaDelta_Ct)
  conds <- gene_df %>% pull(Condition) %>% unique() %>% na.omit()
  conds_to_test <- setdiff(conds, reference_condition)
  
  for (cond in conds_to_test) {
    cond_vals <- gene_df %>% filter(Condition == cond, !Any_ND, !is.na(DeltaDelta_Ct)) %>% pull(DeltaDelta_Ct)
    if (length(ko_vals) >= 2 && length(cond_vals) >= 2) {
      tt <- t.test(cond_vals, ko_vals, var.equal = TRUE, alternative = "two.sided")
      tt_tidy <- broom::tidy(tt) %>% mutate(Gene = g, Condition = cond)
    } else {
      tt_tidy <- tibble(
        estimate = NA_real_, estimate1 = NA_real_, estimate2 = NA_real_,
        statistic = NA_real_, p.value = NA_real_, parameter = NA_real_,
        conf.low = NA_real_, conf.high = NA_real_,
        method = "Two Sample t-test (equal variance)", alternative = "two.sided",
        Gene = g, Condition = cond
      )
    }
    pairwise_list[[length(pairwise_list) + 1]] <- tt_tidy
  }
}

anova_results <- bind_rows(anova_list) %>% select(Gene, everything())
pairwise_vs_KO_ttests <- bind_rows(pairwise_list) %>% relocate(Gene, Condition)

# BH adjustment across all pairwise tests
pairwise_vs_KO_ttests <- pairwise_vs_KO_ttests %>%
  mutate(p_value_adj_BH = ifelse(is.na(p.value), NA_real_, p.adjust(p.value, method = "BH")))

# --------------- WRITE OUTPUT ---------------
sheets <- list()
for (g in names(results_list)) {
  sheet_name <- paste0(substr(g, 1, 25), "_RQ")
  sheets[[sheet_name]] <- results_list[[g]]
}
sheets[["combined"]] <- combined
sheets[["summary"]]  <- summary_df
sheets[["anova_results"]] <- anova_results
sheets[["pairwise_vs_KO_ttests"]] <- pairwise_vs_KO_ttests
sheets[["ko_check_meanRQ"]] <- ko_check
sheets[["raw_collapsed"]] <- collapsed
sheets[["excluded_samples_outlier"]] <- tibble(Sample = excluded_for_outlier)
sheets[["raw_with_flags"]] <- raw %>% select(Sample, Target, CT_raw, undet, CT_num)

write_xlsx(sheets, out_path)
message("Written Excel to: ", out_path)
