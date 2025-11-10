#!/usr/bin/env Rscript

# Quick diagnostic to check LAVA output files

library(data.table)

base_dir <- "/home/lchang24/projects/def-gsarah/lchang24/github/ldsc-lava-kidney-neuro3/lava-analysis/results/sex_specific"
pheno_pairs <- c("egfr_pd", "hematuria_pd", "uacr_pd")
sexes <- c("females", "males")

cat("LAVA OUTPUT FILE DIAGNOSTIC\n")
cat("="*70, "\n\n")

for (sex in sexes) {
  for (pair in pheno_pairs) {
    pair_dir <- file.path(base_dir, sex, "pairwise", pair)
    bivar_file <- file.path(pair_dir, paste0(pair, "_", sex, ".bivar.tsv"))
    univ_file <- file.path(pair_dir, paste0(pair, "_", sex, ".univ.tsv"))
    
    cat(paste0(sex, " - ", pair, ":\n"))
    
    # Check bivariate file
    if (file.exists(bivar_file)) {
      size <- file.info(bivar_file)$size
      if (size > 0) {
        dt <- try(fread(bivar_file, nrows = 1), silent = TRUE)
        if (inherits(dt, "try-error")) {
          cat("  bivar: ERROR reading file\n")
        } else {
          n_rows <- nrow(fread(bivar_file))
          cat(paste("  bivar: OK (", n_rows, "rows,", size, "bytes )\n"))
        }
      } else {
        cat("  bivar: EMPTY FILE\n")
      }
    } else {
      cat("  bivar: FILE NOT FOUND\n")
    }
    
    # Check univariate file
    if (file.exists(univ_file)) {
      size <- file.info(univ_file)$size
      if (size > 0) {
        dt <- try(fread(univ_file, nrows = 1), silent = TRUE)
        if (inherits(dt, "try-error")) {
          cat("  univ: ERROR reading file\n")
        } else {
          n_rows <- nrow(fread(univ_file))
          cat(paste("  univ: OK (", n_rows, "rows,", size, "bytes )\n"))
        }
      } else {
        cat("  univ: EMPTY FILE\n")
      }
    } else {
      cat("  univ: FILE NOT FOUND\n")
    }
    
    cat("\n")
  }
}
