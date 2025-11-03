# Description: Run pairwise LAVA analyses
# 1. hematuria vs PD
# 2. eGFR vs PD
# 3. uacr vs PD

# Load packages -----------------------------------------------------------
library(LAVA)
library(tidyverse)
library(data.table)

# Set arguments -----------------------------------------------------------
project_dir <- "/home/lchang24/projects/def-gsarah/lchang24/github/ldsc-lava-kidney-neuro3/lava-analysis"

# Define paths
ref_prefix <- "/home/lchang24/projects/def-gsarah/sumstats/kidney_neurodegen/lava_ld_ref/ukb_eur/lava-ukb-v1.1"
loc_file <- "/home/lchang24/projects/def-gsarah/sumstats/kidney_neurodegen/lava_locus_file/blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile"
formatted_sumstats_dir <- file.path(project_dir, "output")
output_base_dir <- file.path(project_dir, "results", "pairwise")

# Create output directory
dir.create(output_base_dir, recursive = TRUE, showWarnings = FALSE)

# Load loci ---------------------------------------------------------------
print("Loading loci file...")
loci <- LAVA::read.loci(loc_file)
n_loci <- nrow(loci)
print(paste("Loaded", n_loci, "loci"))

# Set univariate threshold
univar_threshold <- 0.05 / n_loci
print(paste("Univariate threshold:", univar_threshold))

# Print progress milestones
progress <- quantile(x = 1:n_loci, probs = seq(.05, 1, .05)) %>% ceiling()

# Define trait pairs ------------------------------------------------------
trait_pairs <- list(
  list(name = "hematuria_pd", 
       pheno1 = "hematuria", 
       pheno2 = "pd",
       file1 = "formatted_hematuria_sexcombined.tsv",
       file2 = "formatted_pd_sexcombined.tsv"),
  
  list(name = "egfr_pd", 
       pheno1 = "egfr", 
       pheno2 = "pd",
       file1 = "formatted_egfr_sexcombined.tsv",
       file2 = "formatted_pd_sexcombined.tsv"),
  
  list(name = "uacr_pd", 
       pheno1 = "uacr", 
       pheno2 = "pd",
       file1 = "formatted_uacr_sexcombined.tsv",
       file2 = "formatted_pd_sexcombined.tsv")
)

# Run analysis for each pair ----------------------------------------------
for (pair in trait_pairs) {
  
  cat("\n========================================\n")
  cat(paste("Running LAVA for:", pair$name, "\n"))
  cat(paste("Trait 1:", pair$pheno1, "\n"))
  cat(paste("Trait 2:", pair$pheno2, "\n"))
  cat("========================================\n\n")
  
  # Create output directory for this pair
  output_dir <- file.path(output_base_dir, pair$name)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Create info file -----------------------------------------------------
  info_data <- data.frame(
    phenotype = c(pair$pheno1, pair$pheno2),
    cases = c(NA, NA),
    controls = c(NA, NA),
    filename = c(
      file.path(formatted_sumstats_dir, pair$file1),
      file.path(formatted_sumstats_dir, pair$file2)
    ),
    stringsAsFactors = FALSE
  )
  
  # Check if files exist
  print("Checking if formatted summary statistics exist...")
  all_files_exist <- TRUE
  for (i in 1:nrow(info_data)) {
    if (file.exists(info_data$filename[i])) {
      print(paste("✓", info_data$phenotype[i], "- found"))
    } else {
      print(paste("✗ Error:", info_data$phenotype[i], "not found at", info_data$filename[i]))
      all_files_exist <- FALSE
    }
  }
  
  if (!all_files_exist) {
    cat("\nSkipping", pair$name, "due to missing files\n")
    next
  }
  
  info_file <- file.path(output_dir, paste0("info_", pair$name, ".txt"))
  write.table(info_data, info_file, sep = "\t", row.names = FALSE, quote = FALSE)
  print(paste("Info file created:", info_file))
  
  # Create sample overlap file -------------------------------------------
  # Based on LDSC results showing no overlap
  sample_overlap <- matrix(c(1, 0, 0, 1), nrow = 2, ncol = 2)
  rownames(sample_overlap) <- colnames(sample_overlap) <- c(pair$pheno1, pair$pheno2)
  
  sample_overlap_file <- file.path(output_dir, paste0("sample_overlap_", pair$name, ".txt"))
  write.table(sample_overlap, sample_overlap_file, sep = "\t", quote = FALSE, col.names = NA)
  print(paste("Sample overlap file created:", sample_overlap_file))
  
  # Load input data ------------------------------------------------------
  print(paste("Loading input data with", n_loci, "loci..."))
  input <- LAVA::process.input(
    input.info.file = info_file,
    sample.overlap.file = sample_overlap_file,
    ref.prefix = ref_prefix,
    phenos = c(pair$pheno1, pair$pheno2)
  )
  
  # Main analysis --------------------------------------------------------
  print(paste("Starting LAVA analysis for", n_loci, "loci"))
  
  # Initialize result lists
  univar <- bivar <- list()
  
  # Loop through loci
  for (i in 1:n_loci) {
    
    if (i %in% progress) {
      print(paste0("Progress: ", names(progress[which(progress == i)]), 
                   " (locus ", i, " of ", n_loci, ")"))
    }
    
    # Process locus
    locus <- tryCatch(
      LAVA::process.locus(loci[i, ], input),
      error = function(e) {
        warning(paste("Error processing locus", i, ":", e$message))
        return(NULL)
      }
    )
    
    if (!is.null(locus)) {
      
      # Extract locus info
      loc_info <- data.frame(
        locus = locus$id,
        chr = locus$chr,
        start = locus$start,
        stop = locus$stop,
        n_snps = locus$n.snps,
        n_pcs = locus$K
      )
      
      # Run univariate and bivariate tests
      loc_out <- tryCatch(
        LAVA::run.univ.bivar(locus, univ.thresh = univar_threshold),
        error = function(e) {
          warning(paste("Error running tests for locus", i, ":", e$message))
          return(NULL)
        }
      )
      
      if (!is.null(loc_out)) {
        # Store univariate results
        univar[[i]] <- loc_info %>% dplyr::bind_cols(loc_out$univ)
        
        # Store bivariate results if available
        if (!is.null(loc_out$bivar)) {
          bivar[[i]] <- loc_info %>% dplyr::bind_cols(loc_out$bivar)
        }
      }
    }
  }
  
  # Save results ---------------------------------------------------------
  print(paste("\nSaving results for", pair$name, "..."))
  
  # Combine results
  univar_df <- dplyr::bind_rows(univar)
  bivar_df <- dplyr::bind_rows(bivar)
  
  # Save as RDS
  saveRDS(univar, file = file.path(output_dir, paste0(pair$name, ".univ.rds")))
  saveRDS(bivar, file = file.path(output_dir, paste0(pair$name, ".bivar.rds")))
  
  # Save as TSV
  write.table(univar_df, 
              file = file.path(output_dir, paste0(pair$name, ".univ.tsv")),
              sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(bivar_df,
              file = file.path(output_dir, paste0(pair$name, ".bivar.tsv")),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  # Print summary --------------------------------------------------------
  cat("\n=== SUMMARY:", pair$name, "===\n")
  print(paste("Results saved to:", output_dir))
  print(paste("Number of univariate results:", nrow(univar_df)))
  print(paste("Number of bivariate results:", nrow(bivar_df)))
  
  if (nrow(bivar_df) > 0) {
    print(paste("Significant bivariate results (p < 0.05):", sum(bivar_df$p < 0.05, na.rm = TRUE)))
    print(paste("Significant bivariate results (p < 0.01):", sum(bivar_df$p < 0.01, na.rm = TRUE)))
    print(paste("Significant bivariate results (p < 0.001):", sum(bivar_df$p < 0.001, na.rm = TRUE)))
    
    # Show top results
    if (sum(bivar_df$p < 0.05, na.rm = TRUE) > 0) {
      cat("\nTop 10 significant results:\n")
      bivar_df %>%
        filter(p < 0.05) %>%
        arrange(p) %>%
        head(10) %>%
        select(locus, chr, start, stop, phen1, phen2, rho, p) %>%
        print()
    }
  }
  
  cat("\n")
}

cat("\n========================================\n")
cat("All pairwise analyses complete!\n")
cat("========================================\n")
