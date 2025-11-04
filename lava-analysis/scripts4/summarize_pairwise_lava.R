# Description: Summarize pairwise LAVA results for kidney traits vs PD
# Creates comprehensive summary tables and identifies significant findings
# Handles empty bivariate files gracefully

# Load packages -----------------------------------------------------------
library(tidyverse)
library(data.table)

# Set paths ---------------------------------------------------------------
pairwise_dir <- "/home/lchang24/projects/def-gsarah/lchang24/github/ldsc-lava-kidney-neuro3/lava-analysis/results/pairwise"
output_dir <- file.path(pairwise_dir, "summary")

# Create output directory
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Define trait pairs ------------------------------------------------------
trait_pairs <- c("egfr_pd", "hematuria_pd", "uacr_pd")

# Initialize lists to store results
all_univar <- list()
all_bivar <- list()
sig_bivar <- list()

# Read and process each pair ----------------------------------------------
cat("========================================\n")
cat("PAIRWISE LAVA RESULTS SUMMARY\n")
cat("========================================\n\n")

for (pair in trait_pairs) {
  
  cat(paste0("Processing: ", pair, "\n"))
  cat(paste0(rep("-", 40), collapse = ""), "\n")
  
  pair_dir <- file.path(pairwise_dir, pair)
  
  # Read univariate results
  univar_file <- file.path(pair_dir, paste0(pair, ".univ.tsv"))
  bivar_file <- file.path(pair_dir, paste0(pair, ".bivar.tsv"))
  
  if (file.exists(univar_file)) {
    univar <- fread(univar_file)
    all_univar[[pair]] <- univar %>% mutate(analysis = pair)
    
    cat(paste("  Univariate tests:", nrow(univar), "\n"))
    cat(paste("  Significant (p < 0.05):", sum(univar$p < 0.05, na.rm = TRUE), "\n"))
  } else {
    cat("  WARNING: Univariate file not found\n")
  }
  
  # Read bivariate results with error handling
  if (file.exists(bivar_file)) {
    # Check if file is empty
    file_info <- file.info(bivar_file)
    
    if (file_info$size > 100) {  # File has content (more than just header)
      bivar <- tryCatch({
        fread(bivar_file)
      }, error = function(e) {
        cat("  WARNING: Could not read bivariate file (possibly empty)\n")
        return(NULL)
      })
      
      if (!is.null(bivar) && nrow(bivar) > 0) {
        all_bivar[[pair]] <- bivar %>% mutate(analysis = pair)
        
        n_sig_05 <- sum(bivar$p < 0.05, na.rm = TRUE)
        n_sig_01 <- sum(bivar$p < 0.01, na.rm = TRUE)
        n_sig_001 <- sum(bivar$p < 0.001, na.rm = TRUE)
        
        cat(paste("  Bivariate tests:", nrow(bivar), "\n"))
        cat(paste("  Significant (p < 0.05):", n_sig_05, "\n"))
        cat(paste("  Significant (p < 0.01):", n_sig_01, "\n"))
        cat(paste("  Significant (p < 0.001):", n_sig_001, "\n"))
        
        # Store significant results
        if (n_sig_05 > 0) {
          sig_bivar[[pair]] <- bivar %>% 
            filter(p < 0.05) %>%
            arrange(p) %>%
            mutate(analysis = pair)
          
          cat("\n  Top 5 significant results:\n")
          bivar %>%
            filter(p < 0.05) %>%
            arrange(p) %>%
            head(5) %>%
            select(locus, chr, start, stop, rho, p) %>%
            print()
        } else {
          cat("  No significant bivariate correlations found\n")
        }
      } else {
        cat("  Bivariate tests: 0 (file empty or no results)\n")
        cat("  Note: No loci passed univariate threshold for bivariate testing\n")
      }
    } else {
      cat("  Bivariate tests: 0 (file empty)\n")
      cat("  Note: No loci passed univariate threshold for bivariate testing\n")
    }
  } else {
    cat("  WARNING: Bivariate file not found\n")
  }
  
  cat("\n")
}

# Combine all results -----------------------------------------------------
cat("========================================\n")
cat("COMBINED SUMMARY\n")
cat("========================================\n\n")

# Combine univariate
univar_combined <- bind_rows(all_univar)
cat("UNIVARIATE RESULTS:\n")
cat(paste("  Total tests:", nrow(univar_combined), "\n"))
cat(paste("  Significant (p < 0.05):", sum(univar_combined$p < 0.05, na.rm = TRUE), "\n\n"))

# Summary by analysis and phenotype
cat("By analysis and phenotype:\n")
univar_summary <- univar_combined %>%
  group_by(analysis, phen) %>%
  summarise(
    n_loci = n(),
    n_sig_05 = sum(p < 0.05, na.rm = TRUE),
    n_sig_01 = sum(p < 0.01, na.rm = TRUE),
    min_p = min(p, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(analysis, phen)

print(univar_summary, n = Inf)

# Combine bivariate (only if we have results)
if (length(all_bivar) > 0) {
  bivar_combined <- bind_rows(all_bivar)
  
  cat("\n\nBIVARIATE RESULTS:\n")
  cat(paste("  Total tests:", nrow(bivar_combined), "\n"))
  cat(paste("  Significant (p < 0.05):", sum(bivar_combined$p < 0.05, na.rm = TRUE), "\n"))
  cat(paste("  Significant (p < 0.01):", sum(bivar_combined$p < 0.01, na.rm = TRUE), "\n"))
  cat(paste("  Significant (p < 0.001):", sum(bivar_combined$p < 0.001, na.rm = TRUE), "\n\n"))
  
  # Summary by analysis
  cat("By trait pair:\n")
  bivar_summary <- bivar_combined %>%
    group_by(analysis) %>%
    summarise(
      n_tests = n(),
      n_sig_05 = sum(p < 0.05, na.rm = TRUE),
      n_sig_01 = sum(p < 0.01, na.rm = TRUE),
      n_sig_001 = sum(p < 0.001, na.rm = TRUE),
      min_p = min(p, na.rm = TRUE),
      max_abs_rho = max(abs(rho), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(n_sig_05))
  
  print(bivar_summary)
} else {
  cat("\n\nBIVARIATE RESULTS:\n")
  cat("  No bivariate results available (all files empty)\n")
  
  # Create empty bivar_summary for later use
  bivar_summary <- data.frame(
    analysis = trait_pairs,
    n_tests = 0,
    n_sig_05 = 0,
    n_sig_01 = 0,
    n_sig_001 = 0,
    min_p = NA,
    max_abs_rho = NA
  )
}

# Combine significant results ---------------------------------------------
if (length(sig_bivar) > 0) {
  sig_combined <- bind_rows(sig_bivar)
  
  cat("\n\n========================================\n")
  cat("ALL SIGNIFICANT BIVARIATE RESULTS (p < 0.05)\n")
  cat("========================================\n\n")
  
  cat(paste("Total significant loci:", nrow(sig_combined), "\n\n"))
  
  # Show all significant results sorted by p-value
  sig_combined %>%
    arrange(p) %>%
    select(analysis, locus, chr, start, stop, phen1, phen2, rho, p) %>%
    print(n = Inf)
  
  # Chromosome distribution
  cat("\n\nChromosome distribution of significant loci:\n")
  sig_combined %>%
    group_by(analysis, chr) %>%
    summarise(n_loci = n(), .groups = "drop") %>%
    arrange(analysis, chr) %>%
    print(n = Inf)
  
  # Effect direction summary
  cat("\n\nEffect direction (sign of rho):\n")
  sig_combined %>%
    mutate(direction = ifelse(rho > 0, "Positive", "Negative")) %>%
    group_by(analysis, direction) %>%
    summarise(n_loci = n(), .groups = "drop") %>%
    print()
  
} else {
  cat("\n\nNo significant bivariate results found across all analyses.\n")
}

# Save combined results ---------------------------------------------------
cat("\n\n========================================\n")
cat("SAVING RESULTS\n")
cat("========================================\n\n")

# Save all univariate results
write.table(univar_combined, 
            file.path(output_dir, "all_pairwise_univariate.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat(paste("✓ Saved:", file.path(output_dir, "all_pairwise_univariate.tsv"), "\n"))

# Save all bivariate results (if any)
if (length(all_bivar) > 0) {
  write.table(bivar_combined, 
              file.path(output_dir, "all_pairwise_bivariate.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  cat(paste("✓ Saved:", file.path(output_dir, "all_pairwise_bivariate.tsv"), "\n"))
}

# Save significant bivariate results (if any)
if (length(sig_bivar) > 0) {
  write.table(sig_combined, 
              file.path(output_dir, "significant_pairwise_bivariate.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  cat(paste("✓ Saved:", file.path(output_dir, "significant_pairwise_bivariate.tsv"), "\n"))
}

# Save summary statistics
write.table(bivar_summary, 
            file.path(output_dir, "pairwise_summary_stats.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat(paste("✓ Saved:", file.path(output_dir, "pairwise_summary_stats.tsv"), "\n"))

# Save univariate summary
write.table(univar_summary, 
            file.path(output_dir, "pairwise_univariate_summary.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat(paste("✓ Saved:", file.path(output_dir, "pairwise_univariate_summary.tsv"), "\n"))

# Create detailed report --------------------------------------------------
cat("\n\n========================================\n")
cat("GENERATING DETAILED REPORT\n")
cat("========================================\n\n")

report_file <- file.path(output_dir, "pairwise_lava_report.txt")
sink(report_file)

cat(paste0(rep("=", 70), collapse = ""), "\n")
cat("PAIRWISE LAVA ANALYSIS REPORT\n")
cat("Kidney Traits vs Parkinson's Disease Neuropathy\n")
cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat(paste0(rep("=", 70), collapse = ""), "\n\n")

cat("ANALYSES PERFORMED:\n")
cat("1. eGFR vs PD\n")
cat("2. Hematuria vs PD\n")
cat("3. uACR vs PD\n\n")

cat(paste0(rep("=", 70), collapse = ""), "\n")
cat("OVERALL SUMMARY\n")
cat(paste0(rep("=", 70), collapse = ""), "\n\n")

cat("UNIVARIATE RESULTS:\n")
cat(sprintf("  Total loci tested: ~%d per trait\n", nrow(univar_combined)/6))
cat(sprintf("  Total tests: %d\n", nrow(univar_combined)))
cat(sprintf("  Significant (p < 0.05): %d\n", sum(univar_combined$p < 0.05, na.rm = TRUE)))
cat(sprintf("  Significance rate: %.2f%%\n\n", 
            100 * sum(univar_combined$p < 0.05, na.rm = TRUE) / nrow(univar_combined)))

if (length(all_bivar) > 0) {
  cat("BIVARIATE RESULTS:\n")
  cat(sprintf("  Total tests: %d\n", nrow(bivar_combined)))
  cat(sprintf("  Significant (p < 0.05): %d\n", sum(bivar_combined$p < 0.05, na.rm = TRUE)))
  cat(sprintf("  Significant (p < 0.01): %d\n", sum(bivar_combined$p < 0.01, na.rm = TRUE)))
  cat(sprintf("  Significant (p < 0.001): %d\n", sum(bivar_combined$p < 0.001, na.rm = TRUE)))
  cat(sprintf("  Significance rate: %.2f%%\n\n", 
              100 * sum(bivar_combined$p < 0.05, na.rm = TRUE) / nrow(bivar_combined)))
} else {
  cat("BIVARIATE RESULTS:\n")
  cat("  No bivariate tests performed (no loci passed univariate threshold)\n\n")
}

cat(paste0(rep("=", 70), collapse = ""), "\n")
cat("RESULTS BY TRAIT PAIR\n")
cat(paste0(rep("=", 70), collapse = ""), "\n\n")

for (pair in trait_pairs) {
  cat(toupper(gsub("_", " vs ", pair)), "\n")
  cat(paste0(rep("-", 70), collapse = ""), "\n")
  
  if (length(all_bivar) > 0 && pair %in% names(all_bivar)) {
    pair_bivar <- bivar_combined %>% filter(analysis == pair)
    n_sig <- sum(pair_bivar$p < 0.05, na.rm = TRUE)
    
    cat(sprintf("  Loci tested: %d\n", nrow(pair_bivar)))
    cat(sprintf("  Significant: %d (p < 0.05)\n", n_sig))
    
    if (n_sig > 0) {
      sig_loci <- pair_bivar %>% filter(p < 0.05) %>% arrange(p)
      
      cat(sprintf("  Most significant: p = %.2e\n", min(sig_loci$p)))
      cat(sprintf("  Strongest correlation: rho = %.3f\n", 
                  sig_loci$rho[which.max(abs(sig_loci$rho))]))
      
      cat("\n  Top significant loci:\n")
      sig_loci %>%
        head(10) %>%
        mutate(
          position = paste0("chr", chr, ":", start, "-", stop),
          rho_ci = sprintf("%.3f [%.3f, %.3f]", rho, rho.lower, rho.upper)
        ) %>%
        select(locus, position, rho_ci, p) %>%
        print()
    } else {
      cat("  No significant correlations found.\n")
    }
  } else {
    cat("  No bivariate tests performed (no loci passed univariate threshold).\n")
  }
  
  cat("\n")
}

if (length(sig_bivar) > 0) {
  cat(paste0(rep("=", 70), collapse = ""), "\n")
  cat("KEY FINDINGS\n")
  cat(paste0(rep("=", 70), collapse = ""), "\n\n")
  
  # Identify patterns
  chr_distribution <- sig_combined %>%
    group_by(chr) %>%
    summarise(n = n(), .groups = "drop") %>%
    arrange(desc(n))
  
  if (nrow(chr_distribution) > 0) {
    cat("CHROMOSOMAL CLUSTERING:\n")
    if (chr_distribution$n[1] >= 5) {
      cat(sprintf("  Strong clustering on chromosome %d (%d loci, %.1f%% of significant)\n",
                  chr_distribution$chr[1], 
                  chr_distribution$n[1],
                  100 * chr_distribution$n[1] / nrow(sig_combined)))
    } else {
      cat(sprintf("  Chromosome %d has the most significant loci (%d loci)\n",
                  chr_distribution$chr[1], chr_distribution$n[1]))
    }
    cat("\n")
  }
  
  # Check for MHC region
  mhc_loci <- sig_combined %>%
    filter(chr == 6, start >= 25000000, stop <= 35000000)
  
  if (nrow(mhc_loci) > 0) {
    cat("MHC REGION (chr6: 25-35 Mb):\n")
    cat(sprintf("  %d significant loci in MHC region\n", nrow(mhc_loci)))
    cat(sprintf("  Analyses: %s\n", 
                paste(unique(mhc_loci$analysis), collapse = ", ")))
    cat("\n")
  }
  
  # Effect direction
  cat("EFFECT DIRECTIONS:\n")
  sig_combined %>%
    mutate(
      direction = case_when(
        rho > 0.5 ~ "Strong positive",
        rho > 0 ~ "Positive",
        rho < -0.5 ~ "Strong negative",
        rho < 0 ~ "Negative"
      )
    ) %>%
    group_by(analysis, direction) %>%
    summarise(n = n(), .groups = "drop") %>%
    print()
} else {
  cat(paste0(rep("=", 70), collapse = ""), "\n")
  cat("KEY FINDINGS\n")
  cat(paste0(rep("=", 70), collapse = ""), "\n\n")
  cat("No significant bivariate correlations detected.\n")
  cat("Possible reasons:\n")
  cat("  - No loci passed univariate threshold in both traits\n")
  cat("  - True absence of local genetic correlations\n")
  cat("  - Insufficient power to detect weak correlations\n")
}

cat("\n")
cat(paste0(rep("=", 70), collapse = ""), "\n")
cat("END OF REPORT\n")
cat(paste0(rep("=", 70), collapse = ""), "\n")

sink()

cat(paste("✓ Saved detailed report:", report_file, "\n"))

# Create comparison table -------------------------------------------------
cat("\n\nCreating comparison with LDSC results...\n")

ldsc_results <- data.frame(
  pair = c("eGFR - PD", "Hematuria - PD", "uACR - PD"),
  ldsc_rg = c(0.0262, -0.116, -0.0371),
  ldsc_se = c(0.0324, 0.0806, 0.0526),
  ldsc_p = c(0.419, 0.1501, 0.4808),
  stringsAsFactors = FALSE
)

lava_summary <- bivar_summary %>%
  mutate(
    pair = case_when(
      analysis == "egfr_pd" ~ "eGFR - PD",
      analysis == "hematuria_pd" ~ "Hematuria - PD",
      analysis == "uacr_pd" ~ "uACR - PD"
    )
  ) %>%
  select(pair, n_tests, n_sig_05, n_sig_01, min_p, max_abs_rho)

comparison <- ldsc_results %>%
  left_join(lava_summary, by = "pair") %>%
  mutate(
    ldsc_result = sprintf("rg=%.4f (SE=%.4f), p=%.3f", ldsc_rg, ldsc_se, ldsc_p),
    lava_result = ifelse(n_tests > 0,
                         sprintf("%d/%d sig (p<0.05), min_p=%.2e", n_sig_05, n_tests, min_p),
                         "No tests (no loci passed threshold)")
  )

cat("\n\n========================================\n")
cat("LDSC vs LAVA COMPARISON\n")
cat("========================================\n\n")

comparison %>%
  select(pair, ldsc_result, lava_result) %>%
  print()

write.table(comparison, 
            file.path(output_dir, "ldsc_vs_lava_comparison.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat(paste("\n✓ Saved comparison:", file.path(output_dir, "ldsc_vs_lava_comparison.tsv"), "\n"))

cat("\n\n========================================\n")
cat("SUMMARY COMPLETE!\n")
cat("========================================\n")
cat("\nAll outputs saved to:", output_dir, "\n\n")

cat("Files created:\n")
cat("  1. all_pairwise_univariate.tsv\n")
if (length(all_bivar) > 0) {
  cat("  2. all_pairwise_bivariate.tsv\n")
}
if (length(sig_bivar) > 0) {
  cat("  3. significant_pairwise_bivariate.tsv\n")
}
cat("  4. pairwise_summary_stats.tsv\n")
cat("  5. pairwise_univariate_summary.tsv\n")
cat("  6. pairwise_lava_report.txt\n")
cat("  7. ldsc_vs_lava_comparison.tsv\n\n")