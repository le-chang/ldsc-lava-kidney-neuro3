# Description: Compare LAVA results between 1000G and UKB references
# for PD_females vs uacr_female trait pair

# Load packages -----------------------------------------------------------
library(tidyverse)
library(data.table)
library(ggplot2)

# Set arguments -----------------------------------------------------------
project_dir <- "/home/lchang24/projects/def-gsarah/lchang24/github/ldsc-lava-kidney-neuro3/lava-analysis"
results_dir <- file.path(project_dir, "results")
comparison_dir <- file.path(results_dir, "comparison")

# Create comparison directory
dir.create(comparison_dir, recursive = TRUE, showWarnings = FALSE)

# Load results ------------------------------------------------------------
print("Loading 1000G results...")
univar_1000g <- readRDS(file.path(results_dir, "1000g", "PD_females_uacr_female.univ.1000g.rds")) %>%
  bind_rows() %>%
  mutate(reference = "1000G")

bivar_1000g <- readRDS(file.path(results_dir, "1000g", "PD_females_uacr_female.bivar.1000g.rds")) %>%
  bind_rows() %>%
  mutate(reference = "1000G")

print("Loading UKB results...")
univar_ukb <- readRDS(file.path(results_dir, "ukb", "PD_females_uacr_female.univ.ukb.rds")) %>%
  bind_rows() %>%
  mutate(reference = "UKB")

bivar_ukb <- readRDS(file.path(results_dir, "ukb", "PD_females_uacr_female.bivar.ukb.rds")) %>%
  bind_rows() %>%
  mutate(reference = "UKB")

# Combine results ---------------------------------------------------------
print("Combining results...")
univar_combined <- bind_rows(univar_1000g, univar_ukb)
bivar_combined <- bind_rows(bivar_1000g, bivar_ukb)

# Compare univariate results ----------------------------------------------
print("Comparing univariate results...")

# Merge by locus and phenotype
univar_comparison <- univar_1000g %>%
  dplyr::select(locus, phen, rho, rho.se, p, reference) %>%
  dplyr::rename(rho_1000g = rho, rho.se_1000g = rho.se, p_1000g = p) %>%
  dplyr::select(-reference) %>%
  full_join(
    univar_ukb %>%
      dplyr::select(locus, phen, rho, rho.se, p) %>%
      dplyr::rename(rho_ukb = rho, rho.se_ukb = rho.se, p_ukb = p),
    by = c("locus", "phen")
  )

# Calculate correlation
if(nrow(univar_comparison) > 0 & sum(!is.na(univar_comparison$rho_1000g) & !is.na(univar_comparison$rho_ukb)) > 1) {
  cor_univar <- cor(univar_comparison$rho_1000g, univar_comparison$rho_ukb, use = "complete.obs")
  print(paste("Correlation of univariate rho estimates:", round(cor_univar, 3)))
}

# Save univariate comparison
write.table(univar_comparison,
            file = file.path(comparison_dir, "univariate_comparison.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# Compare bivariate results -----------------------------------------------
print("Comparing bivariate results...")

if (nrow(bivar_1000g) > 0 && nrow(bivar_ukb) > 0) {
  # Merge by locus and phenotype pair
  bivar_comparison <- bivar_1000g %>%
    dplyr::select(locus, phen1, phen2, rho, rho.se, p, reference) %>%
    dplyr::rename(rho_1000g = rho, rho.se_1000g = rho.se, p_1000g = p) %>%
    dplyr::select(-reference) %>%
    full_join(
      bivar_ukb %>%
        dplyr::select(locus, phen1, phen2, rho, rho.se, p) %>%
        dplyr::rename(rho_ukb = rho, rho.se_ukb = rho.se, p_ukb = p),
      by = c("locus", "phen1", "phen2")
    )
  
  # Calculate correlation
  if(sum(!is.na(bivar_comparison$rho_1000g) & !is.na(bivar_comparison$rho_ukb)) > 1) {
    cor_bivar <- cor(bivar_comparison$rho_1000g, bivar_comparison$rho_ukb, use = "complete.obs")
    print(paste("Correlation of bivariate rho estimates:", round(cor_bivar, 3)))
  }
  
  # Save bivariate comparison
  write.table(bivar_comparison,
              file = file.path(comparison_dir, "bivariate_comparison.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  # Create scatter plot
  if(sum(!is.na(bivar_comparison$rho_1000g) & !is.na(bivar_comparison$rho_ukb)) > 0) {
    p <- ggplot(bivar_comparison, aes(x = rho_1000g, y = rho_ukb)) +
      geom_point(alpha = 0.6) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
      labs(
        title = "LAVA Bivariate Results: 1000G vs UKB Reference",
        subtitle = "PD_females vs uacr_female",
        x = "Local genetic correlation (1000G)",
        y = "Local genetic correlation (UKB)"
      ) +
      theme_bw() +
      coord_fixed()
    
    ggsave(file.path(comparison_dir, "bivariate_scatter_plot.pdf"), 
           plot = p, width = 8, height = 8)
    print("Scatter plot saved.")
  }
} else {
  print("No bivariate results found for comparison.")
}

# Summary statistics ------------------------------------------------------
print("\n=== SUMMARY ===")
print(paste("1000G - Number of univariate tests:", nrow(univar_1000g)))
print(paste("UKB - Number of univariate tests:", nrow(univar_ukb)))
print(paste("1000G - Number of bivariate tests:", nrow(bivar_1000g)))
print(paste("UKB - Number of bivariate tests:", nrow(bivar_ukb)))

# Significant results (p < 0.05)
if(nrow(bivar_1000g) > 0) {
  print(paste("1000G - Significant bivariate results (p < 0.05):", 
              sum(bivar_1000g$p < 0.05, na.rm = TRUE)))
}
if(nrow(bivar_ukb) > 0) {
  print(paste("UKB - Significant bivariate results (p < 0.05):", 
              sum(bivar_ukb$p < 0.05, na.rm = TRUE)))
}

# Create summary report
summary_report <- data.frame(
  Reference = c("1000G", "UKB"),
  Univariate_Tests = c(nrow(univar_1000g), nrow(univar_ukb)),
  Bivariate_Tests = c(nrow(bivar_1000g), nrow(bivar_ukb)),
  Bivariate_Significant = c(
    sum(bivar_1000g$p < 0.05, na.rm = TRUE),
    sum(bivar_ukb$p < 0.05, na.rm = TRUE)
  )
)

write.table(summary_report,
            file = file.path(comparison_dir, "summary_report.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

print(paste("\nComparison complete! Results saved to", comparison_dir))
