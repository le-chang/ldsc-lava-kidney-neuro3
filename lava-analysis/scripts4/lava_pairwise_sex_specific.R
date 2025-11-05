# Description: Run sex-specific pairwise LAVA analyses
# Separate analyses for females and males
# Pairs: eGFR-PD, hematuria-PD, uACR-PD

# Load packages -----------------------------------------------------------
library(LAVA)
library(tidyverse)
library(data.table)

# Set arguments -----------------------------------------------------------
project_dir <- "/home/lchang24/projects/def-gsarah/lchang24/github/ldsc-lava-kidney-neuro3/lava-analysis"

# Define paths
ref_prefix <- "/home/lchang24/projects/def-gsarah/sumstats/kidney_neurodegen/lava_ld_ref/ukb_eur/lava-ukb-v1.1"
loc_file <- "/home/lchang24/projects/def-gsarah/sumstats/kidney_neurodegen/lava_locus_file/blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile"

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

# Define analyses ---------------------------------------------------------
sexes <- c("females", "males")
trait_pairs <- list(
  list(name = "egfr_pd", kidney = "egfr", neuro = "pd"),
  list(name = "hematuria_pd", kidney = "hematuria", neuro = "pd"),
  list(name = "uacr_pd", kidney = "uacr", neuro = "pd")
)

# Run analyses ------------------------------------------------------------
for (sex in sexes) {
  
  cat("\n========================================\n")
  cat(paste("ANALYZING:", toupper(sex), "\n"))
  cat("========================================\n\n")
  
  formatted_sumstats_dir <- file.path(project_dir, "output", "sex_specific", sex)
  output_base_dir <- file.path(project_dir, "results", "sex_specific", sex, "pairwise")
  
  for (pair in trait_pairs) {
    
    cat(paste0("\nRunning: ", pair$name, " (", sex, ")\n"))
    cat(paste0(rep("-", 60), collapse = ""), "\n")
    
    # Create output directory
    output_dir <- file.path(output_base_dir, pair$name)
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Define file names based on sex
    if (sex == "females") {
      kidney_file <- paste0("formatted_", pair$kidney, "_female.tsv")
      neuro_file <- "formatted_pd_female.tsv"
    } else {
      kidney_file <- paste0("formatted_", pair$kidney, "_male.tsv")
      neuro_file <- "formatted_pd_male.tsv"
    }
    
    # Create info file
    info_data <- data.frame(
      phenotype = c(pair$kidney, pair$neuro),
      cases = c(NA, NA),
      controls = c(NA, NA),
      filename = c(
        file.path(formatted_sumstats_dir, kidney_file),
        file.path(formatted_sumstats_dir, neuro_file)
      ),
      stringsAsFactors = FALSE
    )
    
    # Check if files exist
    print("Checking files...")
    all_files_exist <- TRUE
    for (i in 1:nrow(info_data)) {
      if (file.exists(info_data$filename[i])) {
        print(paste("  ✓", info_data$phenotype[i]))
      } else {
        print(paste("  ✗ Error:", info_data$phenotype[i], "not found"))
        all_files_exist <- FALSE
      }
    }
    
    if (!all_files_exist) {
      cat("Skipping", pair$name, "due to missing files\n")
      next
    }
    
    info_file <- file.path(output_dir, paste0("info_", pair$name, "_", sex, ".txt"))
    write.table(info_data, info_file, sep = "\t", row.names = FALSE, quote = FALSE)
    
    # Create sample overlap file (assume no overlap between kidney and PD)
    sample_overlap <- matrix(c(1, 0, 0, 1), nrow = 2, ncol = 2)
    rownames(sample_overlap) <- colnames(sample_overlap) <- c(pair$kidney, pair$neuro)
    
    sample_overlap_file <- file.path(output_dir, paste0("sample_overlap_", pair$name, "_", sex, ".txt"))
    write.table(sample_overlap, sample_overlap_file, sep = "\t", quote = FALSE, col.names = NA)
    
    # Load input data
    print(paste("Loading input data..."))
    input <- LAVA::process.input(
      input.info.file = info_file,
      sample.overlap.file = sample_overlap_file,
      ref.prefix = ref_prefix,
      phenos = c(pair$kidney, pair$neuro)
    )
    
    # Main analysis
    print(paste("Running LAVA for", n_loci, "loci..."))
    
    univar <- bivar <- list()
    
    for (i in 1:n_loci) {
      
      if (i %in% progress) {
        print(paste0("  Progress: ", names(progress[which(progress == i)]), 
                     " (locus ", i, " of ", n_loci, ")"))
      }
      
      locus <- tryCatch(
        LAVA::process.locus(loci[i, ], input),
        error = function(e) {
          return(NULL)
        }
      )
      
      if (!is.null(locus)) {
        
        loc_info <- data.frame(
          locus = locus$id,
          chr = locus$chr,
          start = locus$start,
          stop = locus$stop,
          n_snps = locus$n.snps,
          n_pcs = locus$K
        )
        
        loc_out <- tryCatch(
          LAVA::run.univ.bivar(locus, univ.thresh = univar_threshold),
          error = function(e) {
            return(NULL)
          }
        )
        
        if (!is.null(loc_out)) {
          univar[[i]] <- loc_info %>% dplyr::bind_cols(loc_out$univ)
          
          if (!is.null(loc_out$bivar)) {
            bivar[[i]] <- loc_info %>% dplyr::bind_cols(loc_out$bivar)
          }
        }
      }
    }
    
    # Save results
    print("Saving results...")
    
    univar_df <- dplyr::bind_rows(univar)
    bivar_df <- dplyr::bind_rows(bivar)
    
    saveRDS(univar, file = file.path(output_dir, paste0(pair$name, "_", sex, ".univ.rds")))
    saveRDS(bivar, file = file.path(output_dir, paste0(pair$name, "_", sex, ".bivar.rds")))
    
    write.table(univar_df, 
                file = file.path(output_dir, paste0(pair$name, "_", sex, ".univ.tsv")),
                sep = "\t", row.names = FALSE, quote = FALSE)
    write.table(bivar_df,
                file = file.path(output_dir, paste0(pair$name, "_", sex, ".bivar.tsv")),
                sep = "\t", row.names = FALSE, quote = FALSE)
    
    # Print summary
    cat("\n=== SUMMARY:", pair$name, "(", sex, ") ===\n")
    print(paste("Univariate tests:", nrow(univar_df)))
    print(paste("Bivariate tests:", nrow(bivar_df)))
    if (nrow(bivar_df) > 0) {
      print(paste("Significant (p<0.05):", sum(bivar_df$p < 0.05, na.rm = TRUE)))
    }
    cat("\n")
  }
}

cat("\n========================================\n")
cat("ALL SEX-SPECIFIC ANALYSES COMPLETE!\n")
cat("========================================\n")
