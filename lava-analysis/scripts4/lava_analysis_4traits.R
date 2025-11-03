# Description: Run LAVA analysis for 4 traits
# uacr, eGFR, hematuria, and PD

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
output_dir <- file.path(project_dir, "results", "4traits")

# Create output directory
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Prepare input info file -------------------------------------------------
info_data <- data.frame(
  phenotype = c("uacr", "egfr", "hematuria", "pd"),
  cases = c(NA, NA, NA, NA),  # NA for quantitative traits (update if binary)
  controls = c(NA, NA, NA, NA),
  filename = c(
    file.path(formatted_sumstats_dir, "formatted_uacr_sexcombined.tsv"),
    file.path(formatted_sumstats_dir, "formatted_egfr_sexcombined.tsv"),
    file.path(formatted_sumstats_dir, "formatted_hematuria_sexcombined.tsv"),
    file.path(formatted_sumstats_dir, "formatted_pd_sexcombined.tsv")
  ),
  stringsAsFactors = FALSE
)

# Check if files exist
print("Checking if formatted summary statistics exist...")
for (i in 1:nrow(info_data)) {
  if (file.exists(info_data$filename[i])) {
    print(paste("✓", info_data$phenotype[i], "- found"))
  } else {
    stop(paste("✗ Error:", info_data$phenotype[i], "not found at", info_data$filename[i]))
  }
}

info_file <- file.path(project_dir, "scripts4", "info_file_4traits.txt")
write.table(info_data, info_file, sep = "\t", row.names = FALSE, quote = FALSE)
print(paste("Info file created:", info_file))

# Create sample overlap file ----------------------------------------------
# 4x4 matrix for sample overlap between traits
# Diagonal = 1 (perfect overlap with self)
# Off-diagonal = sample overlap proportion (0 if no overlap)
# UPDATE these values if you know the actual sample overlap
sample_overlap <- matrix(c(
  1.0, 0.0, 0.0, 0.0,  # uacr
  0.0, 1.0, 0.0, 0.0,  # egfr
  0.0, 0.0, 1.0, 0.0,  # hematuria
  0.0, 0.0, 0.0, 1.0   # pd
), nrow = 4, ncol = 4, byrow = TRUE)

rownames(sample_overlap) <- colnames(sample_overlap) <- c("uacr", "egfr", "hematuria", "pd")

sample_overlap_file <- file.path(project_dir, "scripts4", "sample_overlap_4traits.txt")
write.table(sample_overlap, sample_overlap_file, sep = "\t", quote = FALSE, col.names = NA)
print(paste("Sample overlap file created:", sample_overlap_file))

# Load data ---------------------------------------------------------------
print("Loading loci file...")
loci <- LAVA::read.loci(loc_file)
n_loci <- nrow(loci)

print(paste("Loading input data with", n_loci, "loci..."))
input <- LAVA::process.input(
  input.info.file = info_file,
  sample.overlap.file = sample_overlap_file,
  ref.prefix = ref_prefix,
  phenos = c("uacr", "egfr", "hematuria", "pd")
)

# Main analysis -----------------------------------------------------------
print(paste("Starting LAVA analysis for", n_loci, "loci"))

# Print progress
progress <- quantile(x = 1:n_loci, probs = seq(.05, 1, .05)) %>% ceiling()

# Set univariate threshold
univar_threshold <- 0.05 / n_loci
print(paste("Univariate threshold:", univar_threshold))

# Initialize result lists
univar <- bivar <- list()

# Loop through loci
for (i in 1:n_loci) {
  
  if (i %in% progress) {
    print(paste0("Progress: ", names(progress[which(progress == i)]), " (locus ", i, " of ", n_loci, ")"))
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

# Save results ------------------------------------------------------------
print("Saving results...")

# Combine results
univar_df <- dplyr::bind_rows(univar)
bivar_df <- dplyr::bind_rows(bivar)

# Save as RDS
saveRDS(univar, file = file.path(output_dir, "4traits.univ.rds"))
saveRDS(bivar, file = file.path(output_dir, "4traits.bivar.rds"))

# Save as TSV
write.table(univar_df, 
            file = file.path(output_dir, "4traits.univ.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(bivar_df,
            file = file.path(output_dir, "4traits.bivar.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# Print summary -----------------------------------------------------------
print(paste("\nAnalysis complete! Results saved to", output_dir))
print(paste("Number of univariate results:", nrow(univar_df)))
print(paste("Number of bivariate results:", nrow(bivar_df)))

if (nrow(bivar_df) > 0) {
  print("\n=== Bivariate Results Summary ===")
  print(paste("Total bivariate tests:", nrow(bivar_df)))
  print(paste("Significant results (p < 0.05):", sum(bivar_df$p < 0.05, na.rm = TRUE)))
  print(paste("Significant results (p < 0.01):", sum(bivar_df$p < 0.01, na.rm = TRUE)))
  
  # Summary by trait pair
  if ("phen1" %in% colnames(bivar_df) & "phen2" %in% colnames(bivar_df)) {
    bivar_df %>%
      mutate(pair = paste(phen1, phen2, sep = " - ")) %>%
      group_by(pair) %>%
      summarise(
        n_tests = n(),
        n_sig_05 = sum(p < 0.05, na.rm = TRUE),
        n_sig_01 = sum(p < 0.01, na.rm = TRUE)
      ) %>%
      print()
  }
}