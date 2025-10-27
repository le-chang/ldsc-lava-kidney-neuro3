# Description: run LAVA with 1000G reference for PD_females vs uacr_female
# Load packages -----------------------------------------------------------

# Set R library path
.libPaths("~/.local/R/4.3.1/")

library(LAVA)
library(tidyverse)
library(stringr)

# Set arguments -----------------------------------------------------------

project_dir <- "/home/lchang24/projects/def-gsarah/lchang24/github/ldsc-lava-kidney-neuro3"

# Define phenotypes - using target argument for efficient pairwise testing
phenotypes <- c("PD_females", "uacr_female")

args <- list(
  ref_prefix = file.path(project_dir, "reference_data", "g1000_eur", "g1000_eur"),
  loc_file = file.path(project_dir, "test_loci", "gwas_filtered.loci"),
  info_file = file.path(project_dir, "info_files", "input.info.txt"),
  sample_overlap_file = file.path(project_dir, "sample_overlap", "sample_overlap.txt"),
  phenotypes = phenotypes,
  target = "uacr_female",  # Set uacr_female as target to only test PD_females vs uacr_female
  output_filename = "PD_females_vs_uacr_female_1000G"
)

print("Arguments:")
print(args)

# Load data ---------------------------------------------------------------

cat("\n=== Loading LAVA input data ===\n")
loci <- LAVA::read.loci(args$loc_file)
n_loci <- nrow(loci)

cat(sprintf("Number of loci to test: %d\n", n_loci))

input <- LAVA::process.input(
  input.info.file = args$info_file,
  sample.overlap.file = args$sample_overlap_file,
  ref.prefix = args$ref_prefix,
  phenos = args$phenotypes
)

# Main --------------------------------------------------------------------

# Print progress
cat(sprintf("\n=== Starting LAVA analysis with 1000G reference for %d loci ===\n", n_loci))
progress <- quantile(x = 1:n_loci, probs = seq(.05, 1, .05)) %>% ceiling()

# Set univariate threshold to 0.05/n_loci
univar_threshold <- 0.05 / n_loci
cat(sprintf("Univariate threshold: %.2e\n\n", univar_threshold))

univar <- bivar <- list()

for (i in 1:n_loci) {
  
  if (i %in% progress) {
    cat(sprintf("Progress: %s complete\n", names(progress[which(progress == i)])))
  }
  
  # Process locus
  locus <- tryCatch(
    {
      LAVA::process.locus(loci[i, ], input)
    },
    error = function(e) {
      cat(sprintf("Warning: Error processing locus %d: %s\n", i, e$message))
      return(NULL)
    }
  )
  
  # It is possible that the locus cannot be defined for various reasons (e.g. too few SNPs)
  if (!is.null(locus)) {
    
    # extract some general locus info for the output
    loc_info <- data.frame(
      locus = locus$id,
      chr = locus$chr,
      start = locus$start,
      stop = locus$stop,
      n_snps = locus$n.snps,
      n_pcs = locus$K
    )
    
    # Run the univariate and bivariate tests
    # Using target argument to only test PD_females vs uacr_female
    loc_out <- LAVA::run.univ.bivar(
      locus,
      univ.thresh = univar_threshold,
      target = args$target
    )
    
    # Bind univariate results
    univar[[i]] <- loc_info %>%
      dplyr::bind_cols(loc_out$univ)
    
    # Bind bivariate results if they exist
    if (!is.null(loc_out$bivar)) {
      bivar[[i]] <- loc_info %>%
        dplyr::bind_cols(loc_out$bivar)
    }
    
  }
  
}

# Save data ---------------------------------------------------------------

out_dir <- file.path(project_dir, "sumstats", "kidney_neurodegen", "lava_results", "1000G_ref")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Save as RDS files
saveRDS(
  univar,
  file = file.path(out_dir, paste0(args$output_filename, ".univ.lava.rds"))
)

saveRDS(
  bivar,
  file = file.path(out_dir, paste0(args$output_filename, ".bivar.lava.rds"))
)

# Also save as text files for easier inspection
univar_df <- dplyr::bind_rows(univar)
write.table(
  univar_df,
  file = file.path(out_dir, paste0(args$output_filename, ".univ.lava.txt")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

bivar_df <- dplyr::bind_rows(bivar)
write.table(
  bivar_df,
  file = file.path(out_dir, paste0(args$output_filename, ".bivar.lava.txt")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat(sprintf("\n=== Analysis complete! ===\n"))
cat(sprintf("Results saved to: %s\n", out_dir))
cat(sprintf("Univariate tests: %d loci\n", length(univar[!sapply(univar, is.null)])))
cat(sprintf("Bivariate tests: %d loci\n", length(bivar[!sapply(bivar, is.null)])))
