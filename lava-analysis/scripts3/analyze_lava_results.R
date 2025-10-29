#!/usr/bin/env Rscript
# Description: Robust LAVA results comparison (handles empty files and variable column names)

# Load packages -----------------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(ggplot2)
})

# Set paths ---------------------------------------------------------------
project_dir <- "/home/lchang24/projects/def-gsarah/lchang24/github/ldsc-lava-kidney-neuro3/lava-analysis"
results_dir <- file.path(project_dir, "results")
output_dir <- file.path(results_dir, "comparison")

# Create output directory
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("========================================\n")
cat("LAVA Results Analysis\n")
cat("========================================\n\n")

# Helper function to safely load files ------------------------------------
safe_load <- function(filepath, description) {
  tryCatch({
    if (!file.exists(filepath)) {
      cat(paste("  ", description, "- file not found\n"))
      return(data.frame())
    }
    
    # Check file size
    file_info <- file.info(filepath)
    if (file_info$size < 10) {  # Less than 10 bytes = essentially empty
      cat(paste("  ", description, "- file is empty\n"))
      return(data.frame())
    }
    
    data <- fread(filepath)
    if (nrow(data) == 0) {
      cat(paste("  ", description, "- no data rows\n"))
      return(data.frame())
    }
    
    cat(paste("  ", description, "- loaded", nrow(data), "rows\n"))
    cat(paste("     Columns:", paste(names(data), collapse = ", "), "\n"))
    return(data)
    
  }, error = function(e) {
    cat(paste("  ", description, "- error:", e$message, "\n"))
    return(data.frame())
  })
}

# Load data ---------------------------------------------------------------
cat("Loading data...\n")

univ_1000g <- safe_load(
  file.path(results_dir, "1000g", "PD_females_uacr_female.univ.1000g.tsv"),
  "1000G univariate"
)

univ_ukb <- safe_load(
  file.path(results_dir, "ukb", "PD_females_uacr_female.univ.ukb.tsv"),
  "UKB univariate"
)

bivar_1000g <- safe_load(
  file.path(results_dir, "1000g", "PD_females_uacr_female.bivar.1000g.tsv"),
  "1000G bivariate"
)

bivar_ukb <- safe_load(
  file.path(results_dir, "ukb", "PD_females_uacr_female.bivar.ukb.tsv"),
  "UKB bivariate"
)

cat("\n")

# Check for required columns and standardize names ------------------------
standardize_univ_columns <- function(data, ref_name) {
  if (nrow(data) == 0) return(data)
  
  # Rename h2.obs to h2 for consistency
  if ("h2.obs" %in% names(data)) {
    names(data)[names(data) == "h2.obs"] <- "h2"
  }
  
  return(data)
}

standardize_bivar_columns <- function(data, ref_name) {
  if (nrow(data) == 0) return(data)
  
  # Bivariate columns are already correctly named (rho, p)
  # No changes needed
  return(data)
}

cat("Standardizing column names...\n")
univ_1000g <- standardize_univ_columns(univ_1000g, "1000G univariate")
univ_ukb <- standardize_univ_columns(univ_ukb, "UKB univariate")
bivar_1000g <- standardize_bivar_columns(bivar_1000g, "1000G bivariate")
bivar_ukb <- standardize_bivar_columns(bivar_ukb, "UKB bivariate")
cat("\n")

# Summary statistics ------------------------------------------------------
cat("========================================\n")
cat("SUMMARY STATISTICS\n")
cat("========================================\n\n")

univar_threshold <- 0.05 / 2495
cat(paste("Univariate significance threshold:", format(univar_threshold, scientific = TRUE), "\n"))
cat(paste("(Bonferroni correction: 0.05 / 2495 loci)\n\n"))

# Function to summarize results
summarize_results <- function(univ_data, bivar_data, ref_name) {
  cat(paste(ref_name, "Reference:\n"))
  cat(paste(rep("-", nchar(ref_name) + 11), collapse = ""), "\n")
  
  if (nrow(univ_data) > 0) {
    total_univ <- nrow(univ_data)
    sig_univ <- sum(univ_data$p < univar_threshold, na.rm = TRUE)
    
    cat(paste("  Total univariate tests:", total_univ, "\n"))
    cat(paste("  Significant tests:", sig_univ, 
              paste0("(", round(100*sig_univ/total_univ, 1), "%)")), "\n")
    
    # By phenotype
    for (pheno in c("PD_females", "uacr_female")) {
      pheno_data <- univ_data %>% filter(phen == pheno)
      if (nrow(pheno_data) > 0) {
        sig_pheno <- sum(pheno_data$p < univar_threshold, na.rm = TRUE)
        cat(paste("    ", pheno, ":", nrow(pheno_data), "tests,", 
                  sig_pheno, "significant\n"))
      }
    }
  } else {
    cat("  No univariate data\n")
  }
  
  if (nrow(bivar_data) > 0) {
    sig_bivar <- sum(bivar_data$p < 0.05, na.rm = TRUE)
    cat(paste("  Bivariate tests:", nrow(bivar_data), "\n"))
    cat(paste("  Significant (p<0.05):", sig_bivar, "\n"))
  } else {
    cat("  No bivariate tests performed\n")
  }
  
  cat("\n")
}

summarize_results(univ_1000g, bivar_1000g, "1000G")
summarize_results(univ_ukb, bivar_ukb, "UKB")

# Save summary report -----------------------------------------------------
summary_report <- data.frame(
  Metric = c(
    "Univariate tests (1000G)",
    "Univariate tests (UKB)",
    "Significant univariate (1000G)",
    "Significant univariate (UKB)",
    "Bivariate tests (1000G)",
    "Bivariate tests (UKB)",
    "Significant bivariate (1000G)",
    "Significant bivariate (UKB)"
  ),
  Value = c(
    nrow(univ_1000g),
    nrow(univ_ukb),
    if(nrow(univ_1000g) > 0) sum(univ_1000g$p < univar_threshold) else 0,
    if(nrow(univ_ukb) > 0) sum(univ_ukb$p < univar_threshold) else 0,
    nrow(bivar_1000g),
    nrow(bivar_ukb),
    if(nrow(bivar_1000g) > 0) sum(bivar_1000g$p < 0.05, na.rm = TRUE) else 0,
    if(nrow(bivar_ukb) > 0) sum(bivar_ukb$p < 0.05, na.rm = TRUE) else 0
  )
)

write.table(summary_report,
            file.path(output_dir, "summary_report.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("\n")
print(summary_report)
cat("\n")

# Compare univariate results ----------------------------------------------
if (nrow(univ_1000g) > 0 && nrow(univ_ukb) > 0) {
  cat("========================================\n")
  cat("UNIVARIATE COMPARISON\n")
  cat("========================================\n\n")
  
  # Check if required columns exist
  required_cols_1000g <- c("locus", "phen", "h2", "p")
  required_cols_ukb <- c("locus", "phen", "h2", "p")
  
  missing_1000g <- setdiff(required_cols_1000g, names(univ_1000g))
  missing_ukb <- setdiff(required_cols_ukb, names(univ_ukb))
  
  if (length(missing_1000g) > 0) {
    cat(paste("WARNING: 1000G data missing columns:", paste(missing_1000g, collapse = ", "), "\n"))
  }
  if (length(missing_ukb) > 0) {
    cat(paste("WARNING: UKB data missing columns:", paste(missing_ukb, collapse = ", "), "\n"))
  }
  
  # Proceed with comparison if key columns exist
  if (all(c("locus", "phen", "p") %in% names(univ_1000g)) && 
      all(c("locus", "phen", "p") %in% names(univ_ukb))) {
    
    # Prepare 1000G data - select and rename in one step
    cols_1000g_select <- intersect(c("locus", "chr", "start", "stop", "phen", "h2", "p"), 
                                   names(univ_1000g))
    
    data_1000g <- univ_1000g %>%
      select(all_of(cols_1000g_select))
    
    # Create new column names for renaming (excluding identifier columns)
    cols_to_rename_1000g <- setdiff(cols_1000g_select, c("locus", "chr", "start", "stop", "phen"))
    new_names_1000g <- paste0(cols_to_rename_1000g, "_1000g")
    
    # Rename using direct column name assignment
    for (i in seq_along(cols_to_rename_1000g)) {
      old_name <- cols_to_rename_1000g[i]
      new_name <- new_names_1000g[i]
      names(data_1000g)[names(data_1000g) == old_name] <- new_name
    }
    
    # Prepare UKB data - select and rename in one step
    cols_ukb_select <- intersect(c("locus", "phen", "h2", "p"), 
                                 names(univ_ukb))
    
    data_ukb <- univ_ukb %>%
      select(all_of(cols_ukb_select))
    
    # Create new column names for renaming (excluding identifier columns)
    cols_to_rename_ukb <- setdiff(cols_ukb_select, c("locus", "phen"))
    new_names_ukb <- paste0(cols_to_rename_ukb, "_ukb")
    
    # Rename using direct column name assignment
    for (i in seq_along(cols_to_rename_ukb)) {
      old_name <- cols_to_rename_ukb[i]
      new_name <- new_names_ukb[i]
      names(data_ukb)[names(data_ukb) == old_name] <- new_name
    }
    
    # Merge
    univ_comparison <- data_1000g %>%
      full_join(data_ukb, by = c("locus", "phen"))
    
    write.table(univ_comparison,
                file.path(output_dir, "univariate_comparison.tsv"),
                sep = "\t", row.names = FALSE, quote = FALSE)
    
    # Calculate correlation if h2 columns exist
    if ("h2_1000g" %in% names(univ_comparison) && "h2_ukb" %in% names(univ_comparison)) {
      valid_pairs <- univ_comparison %>%
        filter(!is.na(h2_1000g) & !is.na(h2_ukb))
      
      if (nrow(valid_pairs) > 1) {
        cor_all <- cor(valid_pairs$h2_1000g, valid_pairs$h2_ukb, use = "complete.obs")
        cat(paste("Correlation of h2 estimates:", round(cor_all, 3), "\n\n"))
      }
    }
    
    # Top significant loci
    cat("Top 10 significant loci:\n")
    top_loci <- univ_comparison %>%
      mutate(min_p = pmin(p_1000g, p_ukb, na.rm = TRUE)) %>%
      filter(min_p < univar_threshold) %>%
      arrange(min_p) %>%
      head(10)
    
    if (nrow(top_loci) > 0) {
      print(top_loci)
      write.table(top_loci,
                  file.path(output_dir, "top_significant_loci.tsv"),
                  sep = "\t", row.names = FALSE, quote = FALSE)
    } else {
      cat("No significant loci found.\n")
    }
    cat("\n")
    
    # Create plots
    cat("Creating plots...\n")
    
    if ("h2_1000g" %in% names(univ_comparison) && "h2_ukb" %in% names(univ_comparison)) {
      valid_pairs <- univ_comparison %>%
        filter(!is.na(h2_1000g) & !is.na(h2_ukb))
      
      if (nrow(valid_pairs) > 10) {
        # Scatter plot
        p1 <- ggplot(valid_pairs, aes(x = h2_1000g, y = h2_ukb, color = phen)) +
          geom_point(alpha = 0.5) +
          geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
          facet_wrap(~phen, scales = "free") +
          labs(
            title = "LAVA Univariate: 1000G vs UKB",
            x = "Local SNP-heritability (1000G)",
            y = "Local SNP-heritability (UKB)"
          ) +
          theme_bw() +
          theme(legend.position = "none")
        
        ggsave(file.path(output_dir, "univariate_scatter.pdf"), p1, width = 10, height = 5)
        
        # P-value comparison
        p2 <- ggplot(valid_pairs, aes(x = -log10(p_1000g), y = -log10(p_ukb), color = phen)) +
          geom_point(alpha = 0.5) +
          geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
          geom_hline(yintercept = -log10(univar_threshold), linetype = "dotted") +
          geom_vline(xintercept = -log10(univar_threshold), linetype = "dotted") +
          labs(
            title = "LAVA P-values: 1000G vs UKB",
            x = "-log10(p) 1000G",
            y = "-log10(p) UKB",
            color = "Phenotype"
          ) +
          theme_bw()
        
        ggsave(file.path(output_dir, "univariate_pvalue_comparison.pdf"), p2, width = 8, height = 6)
        cat("  Scatter plots saved\n")
      }
    }
    
    # Manhattan plots
    for (ref in c("1000g", "ukb")) {
      ref_data <- if (ref == "1000g") univ_1000g else univ_ukb
      
      if (nrow(ref_data) > 0 && all(c("chr", "start", "stop", "p") %in% names(ref_data))) {
        p_manhattan <- ref_data %>%
          mutate(pos_mb = (start + stop) / 2 / 1e6) %>%
          ggplot(aes(x = pos_mb, y = -log10(p), color = as.factor(chr))) +
          geom_point(alpha = 0.6, size = 1) +
          geom_hline(yintercept = -log10(univar_threshold), 
                     linetype = "dashed", color = "red") +
          facet_grid(phen ~ chr, scales = "free_x", space = "free_x") +
          labs(
            title = paste("LAVA Results -", toupper(ref)),
            x = "Chromosome",
            y = "-log10(p)"
          ) +
          theme_bw() +
          theme(
            legend.position = "none",
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank()
          )
        
        ggsave(file.path(output_dir, paste0("manhattan_", ref, ".pdf")), 
               p_manhattan, width = 16, height = 6)
        cat(paste("  Manhattan plot saved:", ref, "\n"))
      }
    }
  } else {
    cat("Cannot perform comparison - missing required columns\n")
  }
}

# Bivariate results note --------------------------------------------------
cat("\n========================================\n")
cat("BIVARIATE RESULTS\n")
cat("========================================\n\n")

if (nrow(bivar_1000g) == 0 && nrow(bivar_ukb) == 0) {
  cat("No bivariate tests were performed.\n\n")
  cat("This means no genomic loci showed significant local heritability\n")
  cat("for BOTH PD_females AND uacr_female simultaneously.\n\n")
  cat("Interpretation:\n")
  cat("  - PD and uACR may have genetic signals in different genomic regions\n")
  cat("  - Limited evidence for shared local genetic architecture\n")
  cat("  - This is a valid and meaningful scientific finding\n")
} else {
  if (nrow(bivar_1000g) > 0) {
    cat("1000G bivariate results:\n")
    cols_to_show <- intersect(c("locus", "chr", "start", "stop", "rho", "p"), names(bivar_1000g))
    print(bivar_1000g %>% select(all_of(cols_to_show)))
  }
  if (nrow(bivar_ukb) > 0) {
    cat("\nUKB bivariate results:\n")
    cols_to_show <- intersect(c("locus", "chr", "start", "stop", "rho", "p"), names(bivar_ukb))
    print(bivar_ukb %>% select(all_of(cols_to_show)))
  }
}

# Final summary -----------------------------------------------------------
cat("\n========================================\n")
cat("ANALYSIS COMPLETE\n")
cat("========================================\n\n")

cat("Output files saved to:", output_dir, "\n\n")

cat("Key files:\n")
cat("  - summary_report.tsv\n")
cat("  - univariate_comparison.tsv\n")
if (nrow(univ_1000g) > 0 && nrow(univ_ukb) > 0) {
  cat("  - top_significant_loci.tsv\n")
  cat("  - univariate_scatter.pdf\n")
  cat("  - univariate_pvalue_comparison.pdf\n")
  cat("  - manhattan_1000g.pdf\n")
  cat("  - manhattan_ukb.pdf\n")
}
cat("\n")

cat("Key findings:\n")
cat(paste("  - 1000G:", sum(univ_1000g$p < univar_threshold, na.rm = TRUE), 
          "significant univariate loci\n"))
cat(paste("  - UKB:", sum(univ_ukb$p < univar_threshold, na.rm = TRUE), 
          "significant univariate loci\n"))
cat(paste("  - Bivariate correlations:", nrow(bivar_1000g) + nrow(bivar_ukb), "\n"))
cat("\n")