# Load packages
library(tidyverse)
library(data.table)

# Set paths
results_dir <- "/home/lchang24/projects/def-gsarah/lchang24/github/ldsc-lava-kidney-neuro3/lava-analysis/results/4traits"

# Read results
univar <- fread(file.path(results_dir, "4traits.univ.tsv"))
bivar <- fread(file.path(results_dir, "4traits.bivar.tsv"))

cat("\n=== UNIVARIATE RESULTS SUMMARY ===\n")
cat(paste("Total loci tested:", nrow(univar)/4, "\n"))
cat(paste("Total univariate tests:", nrow(univar), "\n"))
cat(paste("Significant (p < 0.05):", sum(univar$p < 0.05, na.rm = TRUE), "\n"))

# Summary by phenotype
univar %>%
  group_by(phen) %>%
  summarise(
    n_tests = n(),
    n_sig_05 = sum(p < 0.05, na.rm = TRUE),
    n_sig_01 = sum(p < 0.01, na.rm = TRUE)
  ) %>%
  print()

cat("\n=== BIVARIATE RESULTS SUMMARY ===\n")
cat(paste("Total bivariate tests:", nrow(bivar), "\n"))
cat(paste("Significant (p < 0.05):", sum(bivar$p < 0.05, na.rm = TRUE), "\n"))
cat(paste("Significant (p < 0.01):", sum(bivar$p < 0.01, na.rm = TRUE), "\n"))
cat(paste("Significant (p < 0.001):", sum(bivar$p < 0.001, na.rm = TRUE), "\n"))

cat("\n=== RESULTS BY TRAIT PAIR ===\n")
bivar %>%
  mutate(pair = paste(phen1, phen2, sep = " - ")) %>%
  group_by(pair) %>%
  summarise(
    n_tests = n(),
    n_sig_05 = sum(p < 0.05, na.rm = TRUE),
    n_sig_01 = sum(p < 0.01, na.rm = TRUE),
    min_p = min(p, na.rm = TRUE)
  ) %>%
  arrange(desc(n_sig_05)) %>%
  print()

cat("\n=== TOP 20 MOST SIGNIFICANT BIVARIATE RESULTS ===\n")
bivar %>%
  arrange(p) %>%
  head(20) %>%
  select(locus, chr, start, stop, phen1, phen2, rho, p) %>%
  print()

# Save summary
write.table(bivar %>% filter(p < 0.05) %>% arrange(p), 
            file.path(results_dir, "significant_bivar_results.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("\n=== Significant results saved to: ===\n")
cat(paste(file.path(results_dir, "significant_bivar_results.tsv"), "\n"))
