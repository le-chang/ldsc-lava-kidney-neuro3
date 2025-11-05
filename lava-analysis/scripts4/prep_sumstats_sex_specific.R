# Description: Prepare sex-specific summary statistics for LAVA analysis
# Formats female and male data separately for kidney traits and PD

# Load packages -----------------------------------------------------------
library(tidyverse)
library(data.table)

# Set paths ---------------------------------------------------------------
sumstats_dir <- "/home/lchang24/projects/def-gsarah/sumstats/kidney_neurodegen"
output_dir <- "/home/lchang24/projects/def-gsarah/lchang24/github/ldsc-lava-kidney-neuro3/lava-analysis/output/sex_specific"

# Create output directories
dir.create(file.path(output_dir, "females"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "males"), recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# FEMALES
# ==============================================================================

cat("========================================\n")
cat("FORMATTING FEMALE DATA\n")
cat("========================================\n\n")

# Format PD females -------------------------------------------------------
cat("Processing PD females...\n")
pd_f <- fread(file.path(sumstats_dir, "PDnp_females.lava.gz"))
cat("  Columns:", paste(colnames(pd_f), collapse=", "), "\n")
cat("  Rows:", nrow(pd_f), "\n")

# SNP, CHR, BP, A1, A2, BETA, SE, P, A1_AF, N, variant_ID
pd_f_formatted <- pd_f %>%
  mutate(Z = BETA / SE) %>%
  select(SNP, CHR, BP, A1, A2, N, Z, BETA, SE, P, A1_AF)

fwrite(pd_f_formatted, 
       file.path(output_dir, "females", "formatted_pd_female.tsv"),
       sep = "\t")
cat("  ✓ Saved:", file.path(output_dir, "females", "formatted_pd_female.tsv"), "\n\n")

# Format eGFR females -----------------------------------------------------
cat("Processing eGFR females...\n")
egfr_f <- fread(file.path(sumstats_dir, "eGFRcrea_female_HRC_imputed_regenie_allchrs.regenie"))
cat("  Columns:", paste(colnames(egfr_f), collapse=", "), "\n")
cat("  Rows:", nrow(egfr_f), "\n")

# CHROM, GENPOS, ID, ALLELE0, ALLELE1, A1FREQ, INFO, N, TEST, BETA, SE, CHISQ, LOG10P, EXTRA
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

fwrite(egfr_f_formatted, 
       file.path(output_dir, "females", "formatted_egfr_female.tsv"),
       sep = "\t")
cat("  ✓ Saved:", file.path(output_dir, "females", "formatted_egfr_female.tsv"), "\n\n")

# Format uACR females -----------------------------------------------------
cat("Processing uACR females...\n")
uacr_f <- fread(file.path(sumstats_dir, "uacr_female.tsv"))
cat("  Columns:", paste(colnames(uacr_f), collapse=", "), "\n")
cat("  Rows:", nrow(uacr_f), "\n")

# CHR, BP, SNP, A2, A1, N, BETA, P, MAF, SE
uacr_f_formatted <- uacr_f %>%
  mutate(Z = BETA / SE) %>%
  select(SNP, CHR, BP, A1, A2, N, Z, BETA, SE, P, MAF)

fwrite(uacr_f_formatted, 
       file.path(output_dir, "females", "formatted_uacr_female.tsv"),
       sep = "\t")
cat("  ✓ Saved:", file.path(output_dir, "females", "formatted_uacr_female.tsv"), "\n\n")

# Format hematuria females ------------------------------------------------
cat("Processing hematuria females...\n")
hematuria_f <- fread(file.path(sumstats_dir, "X593.FEMALES.ukb_v3.SAIGE.MAC_20.INFO_0.4.txt.gz"))
cat("  Columns:", paste(colnames(hematuria_f), collapse=", "), "\n")
cat("  Rows:", nrow(hematuria_f), "\n")

# CHR, BP, SNP, other_allele, effect_allele, EAF, BETA, SE, P, N
hematuria_f_formatted <- hematuria_f %>%
  mutate(
    A1 = effect_allele,
    A2 = other_allele,
    Z = BETA / SE,
    MAF = ifelse(EAF > 0.5, 1 - EAF, EAF)
  ) %>%
  select(SNP, CHR, BP, A1, A2, N, Z, BETA, SE, P, MAF)

fwrite(hematuria_f_formatted, 
       file.path(output_dir, "females", "formatted_hematuria_female.tsv"),
       sep = "\t")
cat("  ✓ Saved:", file.path(output_dir, "females", "formatted_hematuria_female.tsv"), "\n\n")

# ==============================================================================
# MALES
# ==============================================================================

cat("========================================\n")
cat("FORMATTING MALE DATA\n")
cat("========================================\n\n")

# Format PD males ---------------------------------------------------------
cat("Processing PD males...\n")
pd_m <- fread(file.path(sumstats_dir, "PDnp_males.lava.gz"))
cat("  Columns:", paste(colnames(pd_m), collapse=", "), "\n")
cat("  Rows:", nrow(pd_m), "\n")

pd_m_formatted <- pd_m %>%
  mutate(Z = BETA / SE) %>%
  select(SNP, CHR, BP, A1, A2, N, Z, BETA, SE, P, A1_AF)

fwrite(pd_m_formatted, 
       file.path(output_dir, "males", "formatted_pd_male.tsv"),
       sep = "\t")
cat("  ✓ Saved:", file.path(output_dir, "males", "formatted_pd_male.tsv"), "\n\n")

# Format eGFR males -------------------------------------------------------
cat("Processing eGFR males...\n")
egfr_m <- fread(file.path(sumstats_dir, "eGFRcrea_male_HRC_imputed_regenie_allchrs.regenie"))
cat("  Columns:", paste(colnames(egfr_m), collapse=", "), "\n")
cat("  Rows:", nrow(egfr_m), "\n")

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

fwrite(egfr_m_formatted, 
       file.path(output_dir, "males", "formatted_egfr_male.tsv"),
       sep = "\t")
cat("  ✓ Saved:", file.path(output_dir, "males", "formatted_egfr_male.tsv"), "\n\n")

# Format uACR males -------------------------------------------------------
cat("Processing uACR males...\n")
uacr_m <- fread(file.path(sumstats_dir, "uacr_male.tsv"))
cat("  Columns:", paste(colnames(uacr_m), collapse=", "), "\n")
cat("  Rows:", nrow(uacr_m), "\n")

uacr_m_formatted <- uacr_m %>%
  mutate(Z = BETA / SE) %>%
  select(SNP, CHR, BP, A1, A2, N, Z, BETA, SE, P, MAF)

fwrite(uacr_m_formatted, 
       file.path(output_dir, "males", "formatted_uacr_male.tsv"),
       sep = "\t")
cat("  ✓ Saved:", file.path(output_dir, "males", "formatted_uacr_male.tsv"), "\n\n")

# Format hematuria males --------------------------------------------------
cat("Processing hematuria males...\n")
hematuria_m <- fread(file.path(sumstats_dir, "X593.MALES.ukb_v3.SAIGE.MAC_20.INFO_0.4.txt.gz"))
cat("  Columns:", paste(colnames(hematuria_m), collapse=", "), "\n")
cat("  Rows:", nrow(hematuria_m), "\n")

hematuria_m_formatted <- hematuria_m %>%
  mutate(
    A1 = effect_allele,
    A2 = other_allele,
    Z = BETA / SE,
    MAF = ifelse(EAF > 0.5, 1 - EAF, EAF)
  ) %>%
  select(SNP, CHR, BP, A1, A2, N, Z, BETA, SE, P, MAF)

fwrite(hematuria_m_formatted, 
       file.path(output_dir, "males", "formatted_hematuria_male.tsv"),
       sep = "\t")
cat("  ✓ Saved:", file.path(output_dir, "males", "formatted_hematuria_male.tsv"), "\n\n")

# Summary -----------------------------------------------------------------
cat("========================================\n")
cat("SUMMARY\n")
cat("========================================\n\n")

cat("FEMALES:\n")
cat("  PD:        N =", unique(pd_f_formatted$N), ", SNPs =", nrow(pd_f_formatted), "\n")
cat("  eGFR:      N =", unique(egfr_f_formatted$N), ", SNPs =", nrow(egfr_f_formatted), "\n")
cat("  uACR:      N =", unique(uacr_f_formatted$N), ", SNPs =", nrow(uacr_f_formatted), "\n")
cat("  Hematuria: N =", unique(hematuria_f_formatted$N), ", SNPs =", nrow(hematuria_f_formatted), "\n\n")

cat("MALES:\n")
cat("  PD:        N =", unique(pd_m_formatted$N), ", SNPs =", nrow(pd_m_formatted), "\n")
cat("  eGFR:      N =", unique(egfr_m_formatted$N), ", SNPs =", nrow(egfr_m_formatted), "\n")
cat("  uACR:      N =", unique(uacr_m_formatted$N), ", SNPs =", nrow(uacr_m_formatted), "\n")
cat("  Hematuria: N =", unique(hematuria_m_formatted$N), ", SNPs =", nrow(hematuria_m_formatted), "\n\n")

cat("✓ All formatting complete!\n")
