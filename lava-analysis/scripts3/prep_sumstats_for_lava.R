# Description: Prepare summary statistics for LAVA analysis
# Formats PD_females.lava.gz and uacr_female.tsv to LAVA-compatible format

# Load packages -----------------------------------------------------------
library(tidyverse)
library(data.table)

# Set paths ---------------------------------------------------------------
sumstats_dir <- "/home/lchang24/projects/def-gsarah/sumstats/kidney_neurodegen"
output_dir <- sumstats_dir

# Format PD_females -------------------------------------------------------
print("Formatting PD_females.lava.gz...")

# Read PD_females data
pd_females <- fread(file.path(sumstats_dir, "PD_females.lava.gz"))

# Check column names
print("PD_females columns:")
print(colnames(pd_females))
print(head(pd_females, 3))

# LAVA expects these columns: SNP, CHR, BP, A1, A2, N, Z (or BETA and SE)
# Your file has: SNP, CHR, BP, A1, A2, BETA, SE, P, A1_AF, N, variant_ID

# Calculate Z-score from BETA and SE
pd_females_formatted <- pd_females %>%
  mutate(Z = BETA / SE) %>%
  dplyr::select(SNP, CHR, BP, A1, A2, N, Z, BETA, SE, P, A1_AF)

# Save formatted file
output_file <- file.path(output_dir, "formatted_PD_females.tsv")
fwrite(pd_females_formatted, output_file, sep = "\t")
print(paste("Formatted PD_females saved to:", output_file))

# Format uacr_female ------------------------------------------------------
print("\nFormatting uacr_female.tsv...")

# Read uacr_female data
uacr_female <- fread("/home/lchang24/projects/def-gsarah/lchang24/github/nf-genetic-correlations_fork/nf-genetic-correlations/data/sumstats_original/uacr_female.tsv")

# Check column names
print("uacr_female columns:")
print(colnames(uacr_female))
print(head(uacr_female, 3))

# Your file has: CHR, BP, SNP, A2, A1, N, BETA, P, MAF, SE
# LAVA expects: SNP, CHR, BP, A1, A2, N, Z (or BETA and SE)

# Calculate Z-score from BETA and SE
uacr_female_formatted <- uacr_female %>%
  mutate(Z = BETA / SE) %>%
  dplyr::select(SNP, CHR, BP, A1, A2, N, Z, BETA, SE, P, MAF)

# Save formatted file
output_file <- file.path(output_dir, "formatted_uacr_female.tsv")
fwrite(uacr_female_formatted, output_file, sep = "\t")
print(paste("Formatted uacr_female saved to:", output_file))

# Summary statistics ------------------------------------------------------
print("\n=== SUMMARY ===")
print(paste("PD_females - Number of SNPs:", nrow(pd_females_formatted)))
print(paste("PD_females - Sample size (N):", unique(pd_females_formatted$N)))
print(paste("uacr_female - Number of SNPs:", nrow(uacr_female_formatted)))
print(paste("uacr_female - Sample size (N):", unique(uacr_female_formatted$N)))

print("\nFormatting complete!")