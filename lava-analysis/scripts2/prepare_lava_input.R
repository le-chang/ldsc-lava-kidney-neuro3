# Description: Prepare input.info.txt and sample_overlap.txt files for LAVA analysis
# This creates the input files needed by LAVA from your summary statistics

# Set R library path
.libPaths("~/.local/R/4.3.1/")

library(tidyverse)

# Set paths
project_dir <- " /home/lchang24/projects/def-gsarah/lchang24/github/ldsc-lava-kidney-neuro3"
sumstats_dir <- file.path("/home/lchang24/projects/def-gsarah/sumstats/kidney_neurodegen")

cat("========================================\n")
cat("Preparing LAVA Input Files\n")
cat("========================================\n\n")

# Define your summary statistics files
# Based on your data:
# PD_females.lava.gz has: SNP, CHR, BP, A1, A2, BETA, SE, P, A1_AF, N, variant_ID
# uacr_female.tsv has: CHR, BP, SNP, A2, A1, N, BETA, P, MAF, SE

sumstats_files <- data.frame(
  phenotype = c("PD_females", "uacr_female"),
  filename = c(
    file.path(sumstats_dir, "PD_females.lava.gz"),
    "/home/lchang24/projects/def-gsarah/lchang24/github/nf-genetic-correlations_fork/nf-genetic-correlations/data/sumstats_original/uacr_female.tsv"
  ),
  stringsAsFactors = FALSE
)

# Create input.info.txt with correct column mappings based on your actual data
input_info <- data.frame(
  phenotype = c("PD_females", "uacr_female"),
  filename = sumstats_files$filename,
  SNP = c("SNP", "SNP"),        # Both have SNP column
  A1 = c("A1", "A1"),            # Both have A1 column  
  A2 = c("A2", "A2"),            # Both have A2 column
  N = c("N", "N"),               # Both have N column
  BETA = c("BETA", "BETA"),      # Both have BETA column
  SE = c("SE", "SE"),            # Both have SE column
  P = c("P", "P"),               # Both have P column
  stringsAsFactors = FALSE
)

# Create output directory
info_dir <- file.path(project_dir, "info_files")
dir.create(info_dir, recursive = TRUE, showWarnings = FALSE)

# Write input.info.txt file
input_file <- file.path(info_dir, "input.info.txt")
write.table(
  input_info,
  file = input_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat(sprintf("✓ Input info file created: %s\n\n", input_file))
print(input_info)

# Create sample_overlap.txt file
# Set overlap to 0 if there's no sample overlap between studies
# Adjust this value if you know there is overlap
sample_overlap <- data.frame(
  phenotype1 = "PD_females",
  phenotype2 = "uacr_female",
  overlap = 0  # Change this if there is known sample overlap
)

sample_overlap_file <- file.path(project_dir, "sample_overlap", "sample_overlap.txt")
dir.create(dirname(sample_overlap_file), recursive = TRUE, showWarnings = FALSE)

write.table(
  sample_overlap,
  file = sample_overlap_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat(sprintf("\n✓ Sample overlap file created: %s\n\n", sample_overlap_file))
print(sample_overlap)

cat("\n========================================\n")
cat("Input file preparation complete!\n")
cat("========================================\n")
cat("\nNote: If there IS sample overlap between PD_females and uacr_female,\n")
cat("please edit the sample_overlap.txt file and set the overlap proportion (0-1).\n")
