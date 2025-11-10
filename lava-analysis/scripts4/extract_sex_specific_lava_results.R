#!/usr/bin/env Rscript

# Extract and Summarize Sex-Specific LAVA Results
# Author: Le Chang
# Date: November 2025

library(tidyverse)
library(data.table)

# Set base directory
base_dir <- "/home/lchang24/projects/def-gsarah/lchang24/github/ldsc-lava-kidney-neuro3/lava-analysis/results/sex_specific"

# Define trait pairs
trait_pairs <- c("egfr_pd", "hematuria_pd", "uacr_pd")
sexes <- c("females", "males")

# Function to read and process univariate results
read_univ_results <- function(file_path, sex, trait_pair) {
  if (file.exists(file_path)) {
    df <- fread(file_path)
    df$sex <- sex
    df$trait_pair <- trait_pair
    return(df)
  } else {
    message(paste("File not found:", file_path))
    return(NULL)
  }
}

# Function to read and process bivariate results
read_bivar_results <- function(file_path, sex, trait_pair) {
  if (file.exists(file_path)) {
    df <- fread(file_path)
    df$sex <- sex
    df$trait_pair <- trait_pair
    return(df)
  } else {
    message(paste("File not found:", file_path))
    return(NULL)
  }
}

# Initialize lists to store results
all_univ_results <- list()
all_bivar_results <- list()

# Loop through sexes and trait pairs
counter <- 1
for (sex in sexes) {
  for (trait_pair in trait_pairs) {
    # Construct file paths
    univ_file <- file.path(base_dir, sex, "pairwise", trait_pair, 
                           paste0(trait_pair, "_", sex, ".univ.tsv"))
    bivar_file <- file.path(base_dir, sex, "pairwise", trait_pair, 
                            paste0(trait_pair, "_", sex, ".bivar.tsv"))
    
    # Read univariate results
    message(paste("Reading:", univ_file))
    univ_data <- read_univ_results(univ_file, sex, trait_pair)
    if (!is.null(univ_data)) {
      all_univ_results[[counter]] <- univ_data
    }
    
    # Read bivariate results
    message(paste("Reading:", bivar_file))
    bivar_data <- read_bivar_results(bivar_file, sex, trait_pair)
    if (!is.null(bivar_data)) {
      all_bivar_results[[counter]] <- bivar_data
    }
    
    counter <- counter + 1
  }
}

# Combine all results
univ_combined <- rbindlist(all_univ_results, fill = TRUE)
bivar_combined <- rbindlist(all_bivar_results, fill = TRUE)

# Save combined results
output_dir <- file.path(base_dir, "combined_results")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

fwrite(univ_combined, file.path(output_dir, "all_univariate_results_sex_specific.tsv"), sep = "\t")
fwrite(bivar_combined, file.path(output_dir, "all_bivariate_results_sex_specific.tsv"), sep = "\t")

message("\n=== UNIVARIATE RESULTS SUMMARY ===")
cat("\nTotal univariate tests:", nrow(univ_combined), "\n")
cat("\nBy sex:\n")
print(table(univ_combined$sex))
cat("\nBy trait pair:\n")
print(table(univ_combined$trait_pair))

if ("p" %in% colnames(univ_combined)) {
  sig_univ <- univ_combined[p < 0.05]
  cat("\nSignificant univariate results (p < 0.05):", nrow(sig_univ), "\n")
  
  if (nrow(sig_univ) > 0) {
    cat("\nBy sex:\n")
    print(table(sig_univ$sex))
    cat("\nBy trait pair:\n")
    print(table(sig_univ$trait_pair))
  }
}

message("\n=== BIVARIATE RESULTS SUMMARY ===")
cat("\nTotal bivariate tests:", nrow(bivar_combined), "\n")
cat("\nBy sex:\n")
print(table(bivar_combined$sex))
cat("\nBy trait pair:\n")
print(table(bivar_combined$trait_pair))

if ("p" %in% colnames(bivar_combined)) {
  sig_bivar <- bivar_combined[p < 0.05]
  cat("\nSignificant bivariate results (p < 0.05):", nrow(sig_bivar), "\n")
  
  if (nrow(sig_bivar) > 0) {
    cat("\nBy sex:\n")
    print(table(sig_bivar$sex))
    cat("\nBy trait pair:\n")
    print(table(sig_bivar$trait_pair))
    
    # Save significant results
    fwrite(sig_bivar, file.path(output_dir, "significant_bivariate_results_sex_specific.tsv"), sep = "\t")
    
    cat("\n=== TOP SIGNIFICANT BIVARIATE RESULTS ===\n")
    if ("rho" %in% colnames(sig_bivar)) {
      top_results <- sig_bivar[order(p)][1:min(20, nrow(sig_bivar))]
      print(top_results[, .(locus, phen1, phen2, rho, rho.lower, rho.upper, p, sex, trait_pair)])
    }
  }
}

message("\n=== Output files saved to: ", output_dir, " ===")
message("Files created:")
message("  - all_univariate_results_sex_specific.tsv")
message("  - all_bivariate_results_sex_specific.tsv")
if (exists("sig_bivar") && nrow(sig_bivar) > 0) {
  message("  - significant_bivariate_results_sex_specific.tsv")
}
