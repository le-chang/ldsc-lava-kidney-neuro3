#!/usr/bin/env Rscript

# LAVA Sex-Specific Results Analysis
# Following Nature Genetics 2022 methodology

library(tidyverse)
library(data.table)

# ============================================================================
# CONFIGURATION
# ============================================================================

# Base directories
base_dir <- "/home/lchang24/projects/def-gsarah/lchang24/github/ldsc-lava-kidney-neuro3/lava-analysis/results/sex_specific"
output_dir <- file.path(dirname(base_dir), "summary_results")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Phenotype pairs
pheno_pairs <- c("egfr_pd", "hematuria_pd", "uacr_pd")

# Sex strata
sexes <- c("females", "males")

# ============================================================================
# FUNCTION: Load and Process Results
# ============================================================================

load_lava_results <- function(sex, pheno_pair, base_dir) {
  
  # Construct paths
  pair_dir <- file.path(base_dir, sex, "pairwise", pheno_pair)
  
  # File names
  bivar_file <- file.path(pair_dir, paste0(pheno_pair, "_", sex, ".bivar.tsv"))
  univ_file <- file.path(pair_dir, paste0(pheno_pair, "_", sex, ".univ.tsv"))
  
  # Check if files exist
  if (!file.exists(bivar_file) | !file.exists(univ_file)) {
    warning(paste("Files not found for", sex, pheno_pair))
    return(NULL)
  }
  
  # Load data
  bivar <- fread(bivar_file)
  univ <- fread(univ_file)
  
  # Add metadata
  bivar$sex <- sex
  bivar$pheno_pair <- pheno_pair
  univ$sex <- sex
  univ$pheno_pair <- pheno_pair
  
  # Parse phenotype names from pheno_pair
  phenos <- strsplit(pheno_pair, "_")[[1]]
  bivar$phen1 <- phenos[1]
  bivar$phen2 <- phenos[2]
  
  return(list(bivar = bivar, univ = univ))
}

# ============================================================================
# STEP 1: Load All Results
# ============================================================================

cat("Loading LAVA results...\n")

all_bivar <- list()
all_univ <- list()

for (sex in sexes) {
  for (pair in pheno_pairs) {
    cat(paste("  Loading:", sex, "-", pair, "\n"))
    
    results <- load_lava_results(sex, pair, base_dir)
    
    if (!is.null(results)) {
      key <- paste(sex, pair, sep = "_")
      all_bivar[[key]] <- results$bivar
      all_univ[[key]] <- results$univ
    }
  }
}

# Combine all results
bivar_combined <- rbindlist(all_bivar, fill = TRUE)
univ_combined <- rbindlist(all_univ, fill = TRUE)

cat(paste("\nTotal bivariate tests:", nrow(bivar_combined), "\n"))
cat(paste("Total univariate tests:", nrow(univ_combined), "\n"))

# ============================================================================
# STEP 2: Calculate Significance Thresholds
# Following Nature Genetics LAVA paper methodology
# ============================================================================

cat("\n" , "="*70, "\n")
cat("CALCULATING SIGNIFICANCE THRESHOLDS\n")
cat("="*70, "\n\n")

# Get number of unique loci tested per phenotype pair and sex
loci_summary <- univ_combined %>%
  group_by(sex, pheno_pair, phen) %>%
  summarise(
    n_loci_tested = n(),
    n_loci_sig_nominal = sum(p < 0.05, na.rm = TRUE),
    .groups = "drop"
  )

# For each phenotype pair and sex, count loci where BOTH phenotypes pass univariate filter
# This follows the LAVA paper: "bivariate analysis was performed only for loci 
# in which both phenotypes exhibited univariate signal at P < 0.05/2,495"

# First, get total number of unique loci across all analyses
total_loci <- length(unique(univ_combined$locus))
cat(paste("Total unique loci tested:", total_loci, "\n"))

# Univariate significance threshold (for filtering)
univ_threshold <- 0.05 / total_loci
cat(paste("Univariate filtering threshold: P <", format(univ_threshold, scientific = TRUE), "\n"))
cat(paste("  (0.05 /", total_loci, "loci )\n\n"))

# Count loci passing univariate filter for each phenotype, sex, and pair
univ_sig <- univ_combined %>%
  filter(p < univ_threshold) %>%
  group_by(sex, pheno_pair, phen) %>%
  summarise(n_sig_loci = n(), .groups = "drop")

# For bivariate tests: count loci where BOTH phenotypes pass univariate filter
loci_both_sig <- univ_combined %>%
  filter(p < univ_threshold) %>%
  group_by(sex, pheno_pair, locus) %>%
  summarise(n_phenos_sig = n(), .groups = "drop") %>%
  filter(n_phenos_sig == 2)  # Both phenotypes significant

bivar_tests_performed <- loci_both_sig %>%
  group_by(sex, pheno_pair) %>%
  summarise(n_bivar_tests = n(), .groups = "drop")

# Total bivariate tests performed (sum across all sex-pheno_pair combinations)
total_bivar_tests <- sum(bivar_tests_performed$n_bivar_tests)
cat(paste("Total bivariate tests performed (after univariate filtering):", total_bivar_tests, "\n"))

# Bonferroni correction for bivariate tests
bivar_threshold <- 0.05 / total_bivar_tests
cat(paste("Bivariate significance threshold: P <", format(bivar_threshold, scientific = TRUE), "\n"))
cat(paste("  (0.05 /", total_bivar_tests, "bivariate tests )\n\n"))

# ============================================================================
# STEP 3: Apply Thresholds and Identify Significant Results
# ============================================================================

cat("="*70, "\n")
cat("FILTERING FOR SIGNIFICANT RESULTS\n")
cat("="*70, "\n\n")

# Univariate significant results
univ_sig_results <- univ_combined %>%
  filter(p < univ_threshold) %>%
  arrange(sex, pheno_pair, p)

cat(paste("Univariate significant loci (P <", format(univ_threshold, scientific = TRUE), "):", 
          nrow(univ_sig_results), "\n"))

# Bivariate significant results
bivar_sig_results <- bivar_combined %>%
  filter(p < bivar_threshold) %>%
  arrange(sex, pheno_pair, p)

cat(paste("Bivariate significant correlations (P <", format(bivar_threshold, scientific = TRUE), "):", 
          nrow(bivar_sig_results), "\n\n"))

# ============================================================================
# STEP 4: Detailed Summary Statistics
# ============================================================================

cat("="*70, "\n")
cat("DETAILED SUMMARY BY SEX AND PHENOTYPE PAIR\n")
cat("="*70, "\n\n")

# Summary by sex and phenotype pair
summary_by_group <- bivar_sig_results %>%
  group_by(sex, pheno_pair, phen1, phen2) %>%
  summarise(
    n_sig_loci = n(),
    n_positive_rg = sum(rho > 0, na.rm = TRUE),
    n_negative_rg = sum(rho < 0, na.rm = TRUE),
    mean_rho = mean(rho, na.rm = TRUE),
    median_rho = median(rho, na.rm = TRUE),
    min_p = min(p, na.rm = TRUE),
    n_r2_CI_includes_1 = sum(r2.lower <= 1 & r2.upper >= 1, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    pct_positive = round(100 * n_positive_rg / n_sig_loci, 1),
    pct_negative = round(100 * n_negative_rg / n_sig_loci, 1),
    pct_CI_includes_1 = round(100 * n_r2_CI_includes_1 / n_sig_loci, 1)
  )

print(summary_by_group)

# ============================================================================
# STEP 5: Sex Comparison
# ============================================================================

cat("\n", "="*70, "\n")
cat("SEX-SPECIFIC COMPARISON\n")
cat("="*70, "\n\n")

# Compare number of significant loci between sexes
sex_comparison <- summary_by_group %>%
  select(sex, pheno_pair, n_sig_loci, mean_rho) %>%
  pivot_wider(
    names_from = sex,
    values_from = c(n_sig_loci, mean_rho),
    names_sep = "_"
  )

print(sex_comparison)

# Identify sex-specific loci (significant in one sex but not the other)
females_loci <- bivar_sig_results %>%
  filter(sex == "females") %>%
  select(locus, pheno_pair) %>%
  mutate(sig_in_females = TRUE)

males_loci <- bivar_sig_results %>%
  filter(sex == "males") %>%
  select(locus, pheno_pair) %>%
  mutate(sig_in_males = TRUE)

sex_specific_loci <- full_join(
  females_loci, males_loci,
  by = c("locus", "pheno_pair")
) %>%
  replace_na(list(sig_in_females = FALSE, sig_in_males = FALSE)) %>%
  mutate(
    sex_specificity = case_when(
      sig_in_females & sig_in_males ~ "Both",
      sig_in_females & !sig_in_males ~ "Female-specific",
      !sig_in_females & sig_in_males ~ "Male-specific",
      TRUE ~ "Neither"
    )
  )

sex_specificity_summary <- sex_specific_loci %>%
  group_by(pheno_pair, sex_specificity) %>%
  summarise(n_loci = n(), .groups = "drop") %>%
  arrange(pheno_pair, sex_specificity)

cat("\nSex-specificity of significant loci:\n")
print(sex_specificity_summary)

# ============================================================================
# STEP 6: Effect Size Analysis
# ============================================================================

cat("\n", "="*70, "\n")
cat("EFFECT SIZE ANALYSIS\n")
cat("="*70, "\n\n")

# Analyze effect sizes (rho) for significant results
effect_size_stats <- bivar_sig_results %>%
  group_by(sex, pheno_pair) %>%
  summarise(
    n = n(),
    mean_abs_rho = mean(abs(rho), na.rm = TRUE),
    median_abs_rho = median(abs(rho), na.rm = TRUE),
    min_rho = min(rho, na.rm = TRUE),
    max_rho = max(rho, na.rm = TRUE),
    q25_rho = quantile(rho, 0.25, na.rm = TRUE),
    q75_rho = quantile(rho, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

print(effect_size_stats)

# ============================================================================
# STEP 7: Save Results
# ============================================================================

cat("\n", "="*70, "\n")
cat("SAVING RESULTS\n")
cat("="*70, "\n\n")

# Save all significant results
fwrite(bivar_sig_results, 
       file.path(output_dir, "bivar_significant_results.tsv"),
       sep = "\t")
cat("Saved: bivar_significant_results.tsv\n")

fwrite(univ_sig_results,
       file.path(output_dir, "univ_significant_results.tsv"),
       sep = "\t")
cat("Saved: univ_significant_results.tsv\n")

# Save summary tables
fwrite(summary_by_group,
       file.path(output_dir, "summary_by_sex_and_pair.tsv"),
       sep = "\t")
cat("Saved: summary_by_sex_and_pair.tsv\n")

fwrite(sex_comparison,
       file.path(output_dir, "sex_comparison_summary.tsv"),
       sep = "\t")
cat("Saved: sex_comparison_summary.tsv\n")

fwrite(sex_specific_loci,
       file.path(output_dir, "sex_specific_loci.tsv"),
       sep = "\t")
cat("Saved: sex_specific_loci.tsv\n")

fwrite(sex_specificity_summary,
       file.path(output_dir, "sex_specificity_summary.tsv"),
       sep = "\t")
cat("Saved: sex_specificity_summary.tsv\n")

fwrite(effect_size_stats,
       file.path(output_dir, "effect_size_statistics.tsv"),
       sep = "\t")
cat("Saved: effect_size_statistics.tsv\n")

# Save threshold information
threshold_info <- data.frame(
  metric = c("Total loci tested",
             "Univariate threshold",
             "Total bivariate tests",
             "Bivariate threshold",
             "Univariate significant",
             "Bivariate significant"),
  value = c(total_loci,
            format(univ_threshold, scientific = TRUE),
            total_bivar_tests,
            format(bivar_threshold, scientific = TRUE),
            nrow(univ_sig_results),
            nrow(bivar_sig_results))
)

fwrite(threshold_info,
       file.path(output_dir, "significance_thresholds.tsv"),
       sep = "\t")
cat("Saved: significance_thresholds.tsv\n")

# ============================================================================
# STEP 8: Generate Manuscript-Ready Summary
# ============================================================================

cat("\n", "="*70, "\n")
cat("GENERATING MANUSCRIPT SUMMARY\n")
cat("="*70, "\n\n")

manuscript_summary <- paste0(
  "LAVA SEX-SPECIFIC ANALYSIS RESULTS SUMMARY\n",
  "="*70, "\n\n",
  "Following the methodology of Werme et al. (Nature Genetics 2022),\n",
  "we performed sex-stratified local genetic correlation analysis.\n\n",
  "SIGNIFICANCE THRESHOLDS:\n",
  "- Total unique loci tested: ", total_loci, "\n",
  "- Univariate filtering threshold: P < ", format(univ_threshold, scientific = TRUE), 
  " (Bonferroni: 0.05/", total_loci, ")\n",
  "- Total bivariate tests (after univariate filtering): ", total_bivar_tests, "\n",
  "- Bivariate significance threshold: P < ", format(bivar_threshold, scientific = TRUE),
  " (Bonferroni: 0.05/", total_bivar_tests, ")\n\n",
  "SIGNIFICANT RESULTS:\n",
  "- Univariate significant loci: ", nrow(univ_sig_results), "\n",
  "- Bivariate significant local correlations: ", nrow(bivar_sig_results), "\n\n",
  "BY SEX:\n"
)

for (sex_name in sexes) {
  n_sig <- sum(bivar_sig_results$sex == sex_name)
  n_loci <- length(unique(bivar_sig_results$locus[bivar_sig_results$sex == sex_name]))
  manuscript_summary <- paste0(
    manuscript_summary,
    "- ", str_to_title(sex_name), ": ", n_sig, 
    " significant correlations across ", n_loci, " loci\n"
  )
}

manuscript_summary <- paste0(
  manuscript_summary,
  "\nBY PHENOTYPE PAIR:\n"
)

for (pair in pheno_pairs) {
  phenos <- strsplit(pair, "_")[[1]]
  for (sex_name in sexes) {
    n_sig <- sum(bivar_sig_results$sex == sex_name & bivar_sig_results$pheno_pair == pair)
    if (n_sig > 0) {
      mean_rho <- mean(bivar_sig_results$rho[bivar_sig_results$sex == sex_name & 
                                               bivar_sig_results$pheno_pair == pair], 
                       na.rm = TRUE)
      manuscript_summary <- paste0(
        manuscript_summary,
        "- ", toupper(phenos[1]), " - ", toupper(phenos[2]), " (", sex_name, "): ",
        n_sig, " loci (mean ρ = ", round(mean_rho, 3), ")\n"
      )
    }
  }
}

# Add sex-specificity information
manuscript_summary <- paste0(
  manuscript_summary,
  "\nSEX-SPECIFIC FINDINGS:\n"
)

for (pair in pheno_pairs) {
  pair_specific <- sex_specific_loci %>% filter(pheno_pair == pair)
  n_female_only <- sum(pair_specific$sex_specificity == "Female-specific")
  n_male_only <- sum(pair_specific$sex_specificity == "Male-specific")
  n_both <- sum(pair_specific$sex_specificity == "Both")
  
  if (n_female_only > 0 | n_male_only > 0 | n_both > 0) {
    phenos <- strsplit(pair, "_")[[1]]
    manuscript_summary <- paste0(
      manuscript_summary,
      "- ", toupper(phenos[1]), " - ", toupper(phenos[2]), ":\n",
      "  * Female-specific: ", n_female_only, " loci\n",
      "  * Male-specific: ", n_male_only, " loci\n",
      "  * Both sexes: ", n_both, " loci\n"
    )
  }
}

cat(manuscript_summary)

writeLines(manuscript_summary, 
           file.path(output_dir, "manuscript_summary.txt"))
cat("\nSaved: manuscript_summary.txt\n")

cat("\n", "="*70, "\n")
cat("ANALYSIS COMPLETE!\n")
cat("All results saved to:", output_dir, "\n")
cat("="*70, "\n")