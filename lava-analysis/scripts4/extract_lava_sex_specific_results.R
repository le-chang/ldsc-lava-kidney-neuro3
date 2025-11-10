#!/usr/bin/env Rscript

# LAVA Sex-Specific Results Analysis
# Following Nature Genetics 2022 methodology
# WITH ERROR HANDLING FOR EMPTY FILES

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
# FUNCTION: Load and Process Results (WITH ERROR HANDLING)
# ============================================================================

load_lava_results <- function(sex, pheno_pair, base_dir) {
  
  # Construct paths
  pair_dir <- file.path(base_dir, sex, "pairwise", pheno_pair)
  
  # File names
  bivar_file <- file.path(pair_dir, paste0(pheno_pair, "_", sex, ".bivar.tsv"))
  univ_file <- file.path(pair_dir, paste0(pheno_pair, "_", sex, ".univ.tsv"))
  
  # Check if files exist
  if (!file.exists(bivar_file)) {
    warning(paste("Bivariate file not found for", sex, pheno_pair))
    return(NULL)
  }
  
  if (!file.exists(univ_file)) {
    warning(paste("Univariate file not found for", sex, pheno_pair))
    return(NULL)
  }
  
  # Check if files are empty
  bivar_size <- file.info(bivar_file)$size
  univ_size <- file.info(univ_file)$size
  
  if (is.na(bivar_size) || bivar_size == 0) {
    warning(paste("Bivariate file is empty for", sex, pheno_pair))
    return(NULL)
  }
  
  if (is.na(univ_size) || univ_size == 0) {
    warning(paste("Univariate file is empty for", sex, pheno_pair))
    return(NULL)
  }
  
  # Try to load data with error handling
  bivar <- tryCatch({
    dt <- fread(bivar_file, showProgress = FALSE)
    if (nrow(dt) == 0) {
      warning(paste("No data in bivariate file for", sex, pheno_pair))
      return(NULL)
    }
    dt
  }, error = function(e) {
    warning(paste("Error reading bivariate file for", sex, pheno_pair, ":", e$message))
    return(NULL)
  })
  
  univ <- tryCatch({
    dt <- fread(univ_file, showProgress = FALSE)
    if (nrow(dt) == 0) {
      warning(paste("No data in univariate file for", sex, pheno_pair))
      return(NULL)
    }
    dt
  }, error = function(e) {
    warning(paste("Error reading univariate file for", sex, pheno_pair, ":", e$message))
    return(NULL)
  })
  
  # Check if loading was successful
  if (is.null(bivar) || is.null(univ)) {
    return(NULL)
  }
  
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
cat("="*70, "\n")

all_bivar <- list()
all_univ <- list()
failed_loads <- data.frame(sex = character(), pheno_pair = character(), reason = character())

for (sex in sexes) {
  for (pair in pheno_pairs) {
    cat(paste("  Loading:", sex, "-", pair, "... "))
    
    results <- load_lava_results(sex, pair, base_dir)
    
    if (!is.null(results)) {
      key <- paste(sex, pair, sep = "_")
      all_bivar[[key]] <- results$bivar
      all_univ[[key]] <- results$univ
      cat("SUCCESS\n")
    } else {
      cat("FAILED\n")
      failed_loads <- rbind(failed_loads, 
                           data.frame(sex = sex, pheno_pair = pair, 
                                    reason = "Empty or unreadable file"))
    }
  }
}

cat("\n")

# Check if we have any data
if (length(all_bivar) == 0) {
  stop("ERROR: No data could be loaded from any files. Please check:\n",
       "  1. File paths are correct\n",
       "  2. LAVA analysis completed successfully\n",
       "  3. Output files contain data\n")
}

# Combine all results
bivar_combined <- rbindlist(all_bivar, fill = TRUE)
univ_combined <- rbindlist(all_univ, fill = TRUE)

cat(paste("Successfully loaded:", length(all_bivar), "phenotype-sex combinations\n"))
cat(paste("Total bivariate tests:", nrow(bivar_combined), "\n"))
cat(paste("Total univariate tests:", nrow(univ_combined), "\n\n"))

# Save information about failed loads
if (nrow(failed_loads) > 0) {
  cat("WARNING: Some files could not be loaded:\n")
  print(failed_loads)
  cat("\n")
  fwrite(failed_loads, file.path(output_dir, "failed_file_loads.tsv"), sep = "\t")
}

# ============================================================================
# STEP 2: Calculate Significance Thresholds
# Following Nature Genetics LAVA paper methodology
# ============================================================================

cat("="*70, "\n")
cat("CALCULATING SIGNIFICANCE THRESHOLDS\n")
cat("="*70, "\n\n")

# Get number of unique loci tested
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

cat("Loci passing univariate threshold by phenotype and sex:\n")
print(as.data.frame(univ_sig))
cat("\n")

# For bivariate tests: count loci where BOTH phenotypes pass univariate filter
loci_both_sig <- univ_combined %>%
  filter(p < univ_threshold) %>%
  group_by(sex, pheno_pair, locus) %>%
  summarise(n_phenos_sig = n(), .groups = "drop") %>%
  filter(n_phenos_sig == 2)  # Both phenotypes significant

bivar_tests_performed <- loci_both_sig %>%
  group_by(sex, pheno_pair) %>%
  summarise(n_bivar_tests = n(), .groups = "drop")

cat("Bivariate tests performed (both phenotypes pass univariate filter):\n")
print(as.data.frame(bivar_tests_performed))
cat("\n")

# Total bivariate tests performed
total_bivar_tests <- sum(bivar_tests_performed$n_bivar_tests)
cat(paste("Total bivariate tests performed:", total_bivar_tests, "\n"))

# Check if we have any bivariate tests
if (total_bivar_tests == 0) {
  stop("ERROR: No loci passed the univariate filtering threshold for both phenotypes.\n",
       "This suggests:\n",
       "  1. The LAVA analysis may not have found significant local heritability\n",
       "  2. Effect sizes may be too small to detect\n",
       "  3. Sample sizes may be insufficient\n")
}

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

# Check if we found any significant results
if (nrow(bivar_sig_results) == 0) {
  cat("WARNING: No bivariate correlations passed the Bonferroni-corrected threshold.\n")
  cat("Consider:\n")
  cat("  1. Using a less stringent threshold (e.g., FDR correction)\n")
  cat("  2. Examining nominally significant results (P < 0.05)\n")
  cat("  3. Increasing sample size or power\n\n")
  
  # Output nominal results for exploration
  bivar_nominal <- bivar_combined %>%
    filter(p < 0.05) %>%
    arrange(sex, pheno_pair, p)
  
  cat(paste("Nominally significant correlations (P < 0.05):", nrow(bivar_nominal), "\n\n"))
  
  if (nrow(bivar_nominal) > 0) {
    fwrite(bivar_nominal, 
           file.path(output_dir, "bivar_nominal_results.tsv"),
           sep = "\t")
    cat("Saved nominally significant results to: bivar_nominal_results.tsv\n")
  }
}

# ============================================================================
# STEP 4: Detailed Summary Statistics (if significant results exist)
# ============================================================================

if (nrow(bivar_sig_results) > 0) {
  
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
  
  print(as.data.frame(summary_by_group))
  
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
  
  print(as.data.frame(sex_comparison))
  
  # Identify sex-specific loci
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
  print(as.data.frame(sex_specificity_summary))
  
  # ============================================================================
  # STEP 6: Effect Size Analysis
  # ============================================================================
  
  cat("\n", "="*70, "\n")
  cat("EFFECT SIZE ANALYSIS\n")
  cat("="*70, "\n\n")
  
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
  
  print(as.data.frame(effect_size_stats))
  
} else {
  summary_by_group <- NULL
  sex_comparison <- NULL
  sex_specific_loci <- NULL
  sex_specificity_summary <- NULL
  effect_size_stats <- NULL
}

# ============================================================================
# STEP 7: Save Results
# ============================================================================

cat("\n", "="*70, "\n")
cat("SAVING RESULTS\n")
cat("="*70, "\n\n")

# Save all significant results (if any)
if (nrow(bivar_sig_results) > 0) {
  fwrite(bivar_sig_results, 
         file.path(output_dir, "bivar_significant_results.tsv"),
         sep = "\t")
  cat("Saved: bivar_significant_results.tsv\n")
  
  if (!is.null(summary_by_group)) {
    fwrite(summary_by_group,
           file.path(output_dir, "summary_by_sex_and_pair.tsv"),
           sep = "\t")
    cat("Saved: summary_by_sex_and_pair.tsv\n")
  }
  
  if (!is.null(sex_comparison)) {
    fwrite(sex_comparison,
           file.path(output_dir, "sex_comparison_summary.tsv"),
           sep = "\t")
    cat("Saved: sex_comparison_summary.tsv\n")
  }
  
  if (!is.null(sex_specific_loci)) {
    fwrite(sex_specific_loci,
           file.path(output_dir, "sex_specific_loci.tsv"),
           sep = "\t")
    cat("Saved: sex_specific_loci.tsv\n")
    
    fwrite(sex_specificity_summary,
           file.path(output_dir, "sex_specificity_summary.tsv"),
           sep = "\t")
    cat("Saved: sex_specificity_summary.tsv\n")
  }
  
  if (!is.null(effect_size_stats)) {
    fwrite(effect_size_stats,
           file.path(output_dir, "effect_size_statistics.tsv"),
           sep = "\t")
    cat("Saved: effect_size_statistics.tsv\n")
  }
}

# Always save univariate results and thresholds
fwrite(univ_sig_results,
       file.path(output_dir, "univ_significant_results.tsv"),
       sep = "\t")
cat("Saved: univ_significant_results.tsv\n")

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

# Save complete datasets for further analysis
fwrite(bivar_combined,
       file.path(output_dir, "all_bivar_results.tsv"),
       sep = "\t")
cat("Saved: all_bivar_results.tsv (all tests, not just significant)\n")

fwrite(univ_combined,
       file.path(output_dir, "all_univ_results.tsv"),
       sep = "\t")
cat("Saved: all_univ_results.tsv (all tests, not just significant)\n")

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
  "DATA LOADING:\n",
  "- Phenotype pairs analyzed: ", length(pheno_pairs), " (", paste(pheno_pairs, collapse = ", "), ")\n",
  "- Sex strata: ", length(sexes), " (", paste(sexes, collapse = ", "), ")\n",
  "- Successfully loaded: ", length(all_bivar), " of ", length(pheno_pairs) * length(sexes), 
  " phenotype-sex combinations\n\n",
  "SIGNIFICANCE THRESHOLDS:\n",
  "- Total unique loci tested: ", total_loci, "\n",
  "- Univariate filtering threshold: P < ", format(univ_threshold, scientific = TRUE), 
  " (Bonferroni: 0.05/", total_loci, ")\n",
  "- Total bivariate tests (after univariate filtering): ", total_bivar_tests, "\n",
  "- Bivariate significance threshold: P < ", format(bivar_threshold, scientific = TRUE),
  " (Bonferroni: 0.05/", total_bivar_tests, ")\n\n",
  "SIGNIFICANT RESULTS:\n",
  "- Univariate significant loci: ", nrow(univ_sig_results), "\n",
  "- Bivariate significant local correlations: ", nrow(bivar_sig_results), "\n\n"
)

if (nrow(bivar_sig_results) > 0) {
  manuscript_summary <- paste0(manuscript_summary, "BY SEX:\n")
  for (sex_name in sexes) {
    n_sig <- sum(bivar_sig_results$sex == sex_name)
    n_loci <- length(unique(bivar_sig_results$locus[bivar_sig_results$sex == sex_name]))
    manuscript_summary <- paste0(
      manuscript_summary,
      "- ", str_to_title(sex_name), ": ", n_sig, 
      " significant correlations across ", n_loci, " loci\n"
    )
  }
  
  manuscript_summary <- paste0(manuscript_summary, "\nBY PHENOTYPE PAIR:\n")
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
  
  if (!is.null(sex_specific_loci)) {
    manuscript_summary <- paste0(manuscript_summary, "\nSEX-SPECIFIC FINDINGS:\n")
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
  }
} else {
  manuscript_summary <- paste0(
    manuscript_summary,
    "\nNOTE: No bivariate correlations passed the Bonferroni-corrected threshold.\n",
    "See bivar_nominal_results.tsv for nominally significant results (P < 0.05).\n"
  )
}

cat(manuscript_summary)

writeLines(manuscript_summary, 
           file.path(output_dir, "manuscript_summary.txt"))
cat("\nSaved: manuscript_summary.txt\n")

cat("\n", "="*70, "\n")
cat("ANALYSIS COMPLETE!\n")
cat("All results saved to:", output_dir, "\n")
cat("="*70, "\n")