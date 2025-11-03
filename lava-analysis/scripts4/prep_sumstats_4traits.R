# Description: Prepare summary statistics for LAVA analysis - 4 traits
# Formats uacr_sex_combined, eGFR, hematuria, and PD for LAVA

# Load packages -----------------------------------------------------------
library(tidyverse)
library(data.table)

# Set paths ---------------------------------------------------------------
sumstats_dir <- "/home/lchang24/projects/def-gsarah/sumstats/kidney_neurodegen"
output_dir <- "/home/lchang24/projects/def-gsarah/lchang24/github/ldsc-lava-kidney-neuro3/lava-analysis/output"

# Create output directory if it doesn't exist
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Format uacr_sex_combined ------------------------------------------------
print("Formatting uacr_sex_combined.tsv...")

uacr <- fread(file.path(sumstats_dir, "uacr_sex_combined.tsv"))
print("uacr_sex_combined columns:")
print(colnames(uacr))
print(head(uacr, 3))

# CHR, BP, SNP, A2, A1, N, BETA, P, MAF, SE
# Calculate Z-score from BETA and SE
uacr_formatted <- uacr %>%
  mutate(Z = BETA / SE) %>%
  dplyr::select(SNP, CHR, BP, A1, A2, N, Z, BETA, SE, P, MAF)

output_file <- file.path(output_dir, "formatted_uacr_sexcombined.tsv")
fwrite(uacr_formatted, output_file, sep = "\t")
print(paste("Formatted uacr_sex_combined saved to:", output_file))

# Format eGFR -------------------------------------------------------------
print("\nFormatting metal_eGFR_meta_ea1...")

egfr <- fread(file.path(sumstats_dir, "metal_eGFR_meta_ea1.TBL.map.annot.gc.gz"))
print("eGFR columns:")
print(colnames(egfr))
print(head(egfr, 3))

# MarkerName, Allele1, Allele2, n, Freq1, Effect, StdErr, P.value, etc.
# Need to parse chr:pos from MarkerName
egfr_formatted <- egfr %>%
  separate(MarkerName, into = c("CHR", "BP", "extra"), sep = ":", remove = FALSE, extra = "merge") %>%
  mutate(
    CHR = as.integer(CHR),
    BP = as.integer(BP),
    SNP = RSID,
    A1 = toupper(Allele1),
    A2 = toupper(Allele2),
    N = n,
    BETA = Effect,
    SE = StdErr.GC,  # Using GC-corrected SE
    P = P.value.GC,  # Using GC-corrected P-value
    Z = BETA / SE,
    MAF = ifelse(Freq1 > 0.5, 1 - Freq1, Freq1)
  ) %>%
  dplyr::select(SNP, CHR, BP, A1, A2, N, Z, BETA, SE, P, MAF) %>%
  filter(!is.na(SNP), !is.na(CHR), !is.na(BP))

output_file <- file.path(output_dir, "formatted_egfr_sexcombined.tsv")
fwrite(egfr_formatted, output_file, sep = "\t")
print(paste("Formatted eGFR saved to:", output_file))

# Format hematuria --------------------------------------------------------
print("\nFormatting hematuria_sexcombined.tsv...")

hematuria <- fread(file.path(sumstats_dir, "processed_sumstats/hematuria_sexcombined.tsv"))
print("hematuria columns:")
print(colnames(hematuria))
print(head(hematuria, 3))

# CHR, BP, SNP, A1, A2, BETA, SE, P, N
# Calculate Z-score from BETA and SE
hematuria_formatted <- hematuria %>%
  mutate(Z = BETA / SE) %>%
  dplyr::select(SNP, CHR, BP, A1, A2, N, Z, BETA, SE, P)

output_file <- file.path(output_dir, "formatted_hematuria_sexcombined.tsv")
fwrite(hematuria_formatted, output_file, sep = "\t")
print(paste("Formatted hematuria saved to:", output_file))

# Format PD ---------------------------------------------------------------
print("\nFormatting PDnp_sexcombined.tsv.gz...")

pd <- fread(file.path(sumstats_dir, "PDnp_sexcombined.tsv.gz"))
print("PD columns:")
print(colnames(pd))
print(head(pd, 3))

# chromosome, base_pair_location, effect_allele, other_allele, beta, standard_error, p_value, rsid, etc.
pd_formatted <- pd %>%
  mutate(
    SNP = rsid,
    CHR = chromosome,
    BP = base_pair_location,
    A1 = effect_allele,
    A2 = other_allele,
    BETA = beta,
    SE = standard_error,
    P = p_value,
    Z = BETA / SE,
    MAF = ifelse(effect_allele_frequency > 0.5, 1 - effect_allele_frequency, effect_allele_frequency)
  ) %>%
  dplyr::select(SNP, CHR, BP, A1, A2, Z, BETA, SE, P, MAF) %>%
  filter(!is.na(SNP), !is.na(CHR), !is.na(BP))

# Note: PD doesn't have N column in the header, will need to add it manually
# Based on typical PD GWAS, adding approximate N
pd_formatted$N <- 1000000  # UPDATE THIS with actual sample size if known

output_file <- file.path(output_dir, "formatted_pd_sexcombined.tsv")
fwrite(pd_formatted, output_file, sep = "\t")
print(paste("Formatted PD saved to:", output_file))

# Summary statistics ------------------------------------------------------
print("\n=== SUMMARY ===")
print(paste("uacr_sex_combined - Number of SNPs:", nrow(uacr_formatted)))
print(paste("uacr_sex_combined - Sample size (N):", unique(uacr_formatted$N)))

print(paste("eGFR - Number of SNPs:", nrow(egfr_formatted)))
print(paste("eGFR - Sample size range:", min(egfr_formatted$N), "-", max(egfr_formatted$N)))

print(paste("hematuria - Number of SNPs:", nrow(hematuria_formatted)))
print(paste("hematuria - Sample size (N):", unique(hematuria_formatted$N)))

print(paste("PD - Number of SNPs:", nrow(pd_formatted)))
print(paste("PD - Sample size (N):", unique(pd_formatted$N)))

print("\nFormatting complete!")
