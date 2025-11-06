# Description: Prepare sex-specific summary statistics for LAVA analysis
# Formats female and male data separately for kidney traits and PD
# MEMORY OPTIMIZED VERSION

# Load packages -----------------------------------------------------------
library(tidyverse)
library(data.table)

# Memory optimization settings --------------------------------------------
options(datatable.verbose = FALSE)
gc()  # Initial garbage collection

# Set paths ---------------------------------------------------------------
sumstats_dir <- "/home/lchang24/projects/def-gsarah/sumstats/kidney_neurodegen"
output_dir <- "/home/lchang24/projects/def-gsarah/lchang24/github/ldsc-lava-kidney-neuro3/lava-analysis/output/sex_specific"

# Create output directories
dir.create(file.path(output_dir, "females"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "males"), recursive = TRUE, showWarnings = FALSE)

# Helper function to report memory usage ----------------------------------
report_memory <- function(label = "") {
  mem_used <- gc()[2, 2]  # Get memory used in MB
  cat(sprintf("  [Memory: %.1f MB used] %s\n", mem_used, label))
}

# ==============================================================================
# FEMALES
# ==============================================================================

cat("========================================\n")
cat("FORMATTING FEMALE DATA\n")
cat("========================================\n\n")

# Format PD females -------------------------------------------------------
cat("Processing PD females...\n")
report_memory("Before loading")

pd_f <- fread(
  file.path(sumstats_dir, "PDnp_females.lava.gz"),
  select = c("SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "P", "A1_AF", "N")
)
cat("  Columns:", paste(colnames(pd_f), collapse=", "), "\n")
cat("  Rows:", nrow(pd_f), "\n")
report_memory("After loading")

pd_f_formatted <- pd_f %>%
  mutate(Z = BETA / SE) %>%
  select(SNP, CHR, BP, A1, A2, N, Z, BETA, SE, P, A1_AF)

# Store stats BEFORE removing
female_pd_n <- unique(pd_f_formatted$N)[1]
female_pd_snps <- nrow(pd_f_formatted)

rm(pd_f)
gc()

fwrite(pd_f_formatted, 
       file.path(output_dir, "females", "formatted_pd_female.tsv"),
       sep = "\t")
cat("  ✓ Saved:", file.path(output_dir, "females", "formatted_pd_female.tsv"), "\n")
report_memory("After saving")

rm(pd_f_formatted)
gc()
cat("\n")

# Format eGFR females -----------------------------------------------------
cat("Processing eGFR females...\n")
report_memory("Before loading")

egfr_f <- fread(
  file.path(sumstats_dir, "eGFRcrea_female_HRC_imputed_regenie_allchrs.regenie"),
  select = c("CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1", "A1FREQ", "N", "BETA", "SE", "LOG10P")
)
cat("  Columns:", paste(colnames(egfr_f), collapse=", "), "\n")
cat("  Rows:", nrow(egfr_f), "\n")
report_memory("After loading")

egfr_f_formatted <- egfr_f %>%
  mutate(
    SNP = ID,
    CHR = CHROM,
    BP = GENPOS,
    A1 = ALLELE1,
    A2 = ALLELE0,
    Z = BETA / SE,
    P = 10^(-LOG10P),
    MAF = ifelse(A1FREQ > 0.5, 1 - A1FREQ, A1FREQ)
  ) %>%
  select(SNP, CHR, BP, A1, A2, N, Z, BETA, SE, P, MAF) %>%
  filter(!is.na(SNP), !is.na(CHR), !is.na(BP))

# Store stats BEFORE removing
female_egfr_n <- unique(egfr_f_formatted$N)[1]
female_egfr_snps <- nrow(egfr_f_formatted)

rm(egfr_f)
gc()

fwrite(egfr_f_formatted, 
       file.path(output_dir, "females", "formatted_egfr_female.tsv"),
       sep = "\t")
cat("  ✓ Saved:", file.path(output_dir, "females", "formatted_egfr_female.tsv"), "\n")
report_memory("After saving")

rm(egfr_f_formatted)
gc()
cat("\n")

# Format uACR females -----------------------------------------------------
cat("Processing uACR females...\n")
report_memory("Before loading")

uacr_f <- fread(
  file.path(sumstats_dir, "uacr_female.tsv"),
  select = c("CHR", "BP", "SNP", "A2", "A1", "N", "BETA", "P", "MAF", "SE")
)
cat("  Columns:", paste(colnames(uacr_f), collapse=", "), "\n")
cat("  Rows:", nrow(uacr_f), "\n")
report_memory("After loading")

uacr_f_formatted <- uacr_f %>%
  mutate(Z = BETA / SE) %>%
  select(SNP, CHR, BP, A1, A2, N, Z, BETA, SE, P, MAF)

# Store stats BEFORE removing
female_uacr_n <- unique(uacr_f_formatted$N)[1]
female_uacr_snps <- nrow(uacr_f_formatted)

rm(uacr_f)
gc()

fwrite(uacr_f_formatted, 
       file.path(output_dir, "females", "formatted_uacr_female.tsv"),
       sep = "\t")
cat("  ✓ Saved:", file.path(output_dir, "females", "formatted_uacr_female.tsv"), "\n")
report_memory("After saving")

rm(uacr_f_formatted)
gc()
cat("\n")

# Format hematuria females ------------------------------------------------
cat("Processing hematuria females...\n")
report_memory("Before loading")

hematuria_f <- fread(
  file.path(sumstats_dir, "X593.FEMALES.ukb_v3.SAIGE.MAC_20.INFO_0.4.txt.gz"),
  select = c("CHR", "BP", "SNP", "other_allele", "effect_allele", "EAF", "BETA", "SE", "P", "N")
)
cat("  Columns:", paste(colnames(hematuria_f), collapse=", "), "\n")
cat("  Rows:", nrow(hematuria_f), "\n")
report_memory("After loading")

hematuria_f_formatted <- hematuria_f %>%
  mutate(
    A1 = effect_allele,
    A2 = other_allele,
    Z = BETA / SE,
    MAF = ifelse(EAF > 0.5, 1 - EAF, EAF)
  ) %>%
  select(SNP, CHR, BP, A1, A2, N, Z, BETA, SE, P, MAF)

# Store stats BEFORE removing
female_hematuria_n <- unique(hematuria_f_formatted$N)[1]
female_hematuria_snps <- nrow(hematuria_f_formatted)

rm(hematuria_f)
gc()

fwrite(hematuria_f_formatted, 
       file.path(output_dir, "females", "formatted_hematuria_female.tsv"),
       sep = "\t")
cat("  ✓ Saved:", file.path(output_dir, "females", "formatted_hematuria_female.tsv"), "\n")
report_memory("After saving")

rm(hematuria_f_formatted)
gc()
cat("\n")

cat("Female data processing complete. Clearing memory before males...\n")
gc()
Sys.sleep(2)  # Brief pause to ensure memory is released

# ==============================================================================
# MALES
# ==============================================================================

cat("========================================\n")
cat("FORMATTING MALE DATA\n")
cat("========================================\n\n")

# Format PD males ---------------------------------------------------------
cat("Processing PD males...\n")
report_memory("Before loading")

pd_m <- fread(
  file.path(sumstats_dir, "PDnp_males.lava.gz"),
  select = c("SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "P", "A1_AF", "N")
)
cat("  Columns:", paste(colnames(pd_m), collapse=", "), "\n")
cat("  Rows:", nrow(pd_m), "\n")
report_memory("After loading")

pd_m_formatted <- pd_m %>%
  mutate(Z = BETA / SE) %>%
  select(SNP, CHR, BP, A1, A2, N, Z, BETA, SE, P, A1_AF)

# Store stats BEFORE removing
male_pd_n <- unique(pd_m_formatted$N)[1]
male_pd_snps <- nrow(pd_m_formatted)

rm(pd_m)
gc()

fwrite(pd_m_formatted, 
       file.path(output_dir, "males", "formatted_pd_male.tsv"),
       sep = "\t")
cat("  ✓ Saved:", file.path(output_dir, "males", "formatted_pd_male.tsv"), "\n")
report_memory("After saving")

rm(pd_m_formatted)
gc()
cat("\n")

# Format eGFR males -------------------------------------------------------
cat("Processing eGFR males...\n")
report_memory("Before loading")

egfr_m <- fread(
  file.path(sumstats_dir, "eGFRcrea_male_HRC_imputed_regenie_allchrs.regenie"),
  select = c("CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1", "A1FREQ", "N", "BETA", "SE", "LOG10P")
)
cat("  Columns:", paste(colnames(egfr_m), collapse=", "), "\n")
cat("  Rows:", nrow(egfr_m), "\n")
report_memory("After loading")

egfr_m_formatted <- egfr_m %>%
  mutate(
    SNP = ID,
    CHR = CHROM,
    BP = GENPOS,
    A1 = ALLELE1,
    A2 = ALLELE0,
    Z = BETA / SE,
    P = 10^(-LOG10P),
    MAF = ifelse(A1FREQ > 0.5, 1 - A1FREQ, A1FREQ)
  ) %>%
  select(SNP, CHR, BP, A1, A2, N, Z, BETA, SE, P, MAF) %>%
  filter(!is.na(SNP), !is.na(CHR), !is.na(BP))

# Store stats BEFORE removing
male_egfr_n <- unique(egfr_m_formatted$N)[1]
male_egfr_snps <- nrow(egfr_m_formatted)

rm(egfr_m)
gc()

fwrite(egfr_m_formatted, 
       file.path(output_dir, "males", "formatted_egfr_male.tsv"),
       sep = "\t")
cat("  ✓ Saved:", file.path(output_dir, "males", "formatted_egfr_male.tsv"), "\n")
report_memory("After saving")

rm(egfr_m_formatted)
gc()
cat("\n")

# Format uACR males -------------------------------------------------------
cat("Processing uACR males...\n")
report_memory("Before loading")

uacr_m <- fread(
  file.path(sumstats_dir, "uacr_male.tsv"),
  select = c("CHR", "BP", "SNP", "A2", "A1", "N", "BETA", "P", "MAF", "SE")
)
cat("  Columns:", paste(colnames(uacr_m), collapse=", "), "\n")
cat("  Rows:", nrow(uacr_m), "\n")
report_memory("After loading")

uacr_m_formatted <- uacr_m %>%
  mutate(Z = BETA / SE) %>%
  select(SNP, CHR, BP, A1, A2, N, Z, BETA, SE, P, MAF)

# Store stats BEFORE removing
male_uacr_n <- unique(uacr_m_formatted$N)[1]
male_uacr_snps <- nrow(uacr_m_formatted)

rm(uacr_m)
gc()

fwrite(uacr_m_formatted, 
       file.path(output_dir, "males", "formatted_uacr_male.tsv"),
       sep = "\t")
cat("  ✓ Saved:", file.path(output_dir, "males", "formatted_uacr_male.tsv"), "\n")
report_memory("After saving")

rm(uacr_m_formatted)
gc()
cat("\n")

# Format hematuria males --------------------------------------------------
cat("Processing hematuria males...\n")
report_memory("Before loading")

# This is where the crash occurred - being extra careful with memory
hematuria_m <- fread(
  file.path(sumstats_dir, "X593.MALES.ukb_v3.SAIGE.MAC_20.INFO_0.4.txt.gz"),
  select = c("CHR", "BP", "SNP", "other_allele", "effect_allele", "EAF", "BETA", "SE", "P", "N")
)
cat("  Columns:", paste(colnames(hematuria_m), collapse=", "), "\n")
cat("  Rows:", nrow(hematuria_m), "\n")
report_memory("After loading")

hematuria_m_formatted <- hematuria_m %>%
  mutate(
    A1 = effect_allele,
    A2 = other_allele,
    Z = BETA / SE,
    MAF = ifelse(EAF > 0.5, 1 - EAF, EAF)
  ) %>%
  select(SNP, CHR, BP, A1, A2, N, Z, BETA, SE, P, MAF)

# Store stats BEFORE removing
male_hematuria_n <- unique(hematuria_m_formatted$N)[1]
male_hematuria_snps <- nrow(hematuria_m_formatted)

# Clear original immediately
rm(hematuria_m)
gc()
report_memory("After formatting")

fwrite(hematuria_m_formatted, 
       file.path(output_dir, "males", "formatted_hematuria_male.tsv"),
       sep = "\t")
cat("  ✓ Saved:", file.path(output_dir, "males", "formatted_hematuria_male.tsv"), "\n")
report_memory("After saving")

rm(hematuria_m_formatted)
gc()
cat("\n")

# Summary -----------------------------------------------------------------
cat("========================================\n")
cat("SUMMARY\n")
cat("========================================\n\n")

cat("FEMALES:\n")
cat("  PD:        N =", female_pd_n, ", SNPs =", female_pd_snps, "\n")
cat("  eGFR:      N =", female_egfr_n, ", SNPs =", female_egfr_snps, "\n")
cat("  uACR:      N =", female_uacr_n, ", SNPs =", female_uacr_snps, "\n")
cat("  Hematuria: N =", female_hematuria_n, ", SNPs =", female_hematuria_snps, "\n\n")

cat("MALES:\n")
cat("  PD:        N =", male_pd_n, ", SNPs =", male_pd_snps, "\n")
cat("  eGFR:      N =", male_egfr_n, ", SNPs =", male_egfr_snps, "\n")
cat("  uACR:      N =", male_uacr_n, ", SNPs =", male_uacr_snps, "\n")
cat("  Hematuria: N =", male_hematuria_n, ", SNPs =", male_hematuria_snps, "\n\n")

cat("✓ All formatting complete!\n")
report_memory("Final")