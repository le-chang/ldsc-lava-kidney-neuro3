# Description: Run LAVA analysis with UKB reference
# for PD_females vs uacr_female trait pair

# Load packages -----------------------------------------------------------
library(LAVA)
library(tidyverse)
library(data.table)

# Set arguments -----------------------------------------------------------
project_dir <- "/home/lchang24/projects/def-gsarah/lchang24/github/ldsc-lava-kidney-neuro3/lava-analysis"

# Define paths
ref_prefix <- "/home/lchang24/projects/def-gsarah/lchang24/github/nf-genetic-correlations_fork/nf-genetic-correlations/data/ld_reference/lava_ld_ref/ukb_eur/lava-ukb-v1.1"
loc_file <- "/home/lchang24/projects/def-gsarah/lchang24/github/nf-genetic-correlations_fork/nf-genetic-correlations/data/lava_locus_file/blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile"
formatted_sumstats_dir <- file.path(project_dir, "output")  # Updated to output directory
output_dir <- file.path(project_dir, "results", "ukb")

# Create output directory
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Prepare input info file -------------------------------------------------
# Create info file for LAVA pointing to formatted sumstats in output directory
info_data <- data.frame(
  phenotype = c("PD_females", "uacr_female"),
  cases = c(NA, NA),  # NA for quantitative traits
  controls = c(NA, NA),  # NA for quantitative traits
  filename = c(
    file.path(formatted_sumstats_dir, "formatted_PD_females.tsv"),
    file.path(formatted_sumstats_dir, "formatted_uacr_female.tsv")
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

info_file <- file.path(project_dir, "scripts3", "info_file_ukb.txt")
write.table(info_data, info_file, sep = "\t", row.names = FALSE, quote = FALSE)
print(paste("Info file created:", info_file))

# Create sample overlap file ----------------------------------------------
# Assuming no sample overlap between PD_females and uacr_female
# If there is overlap, update the off-diagonal values
sample_overlap <- matrix(c(1, 0, 0, 1), nrow = 2, ncol = 2)
rownames(sample_overlap) <- colnames(sample_overlap) <- c("PD_females", "uacr_female")

sample_overlap_file <- file.path(project_dir, "scripts3", "sample_overlap_ukb.txt")
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
  phenos = c("PD_females", "uacr_female")
)

# Main analysis -----------------------------------------------------------
print(paste("Starting LAVA analysis with UKB reference for", n_loci, "loci"))

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
saveRDS(univar, file = file.path(output_dir, "PD_females_uacr_female.univ.ukb.rds"))
saveRDS(bivar, file = file.path(output_dir, "PD_females_uacr_female.bivar.ukb.rds"))

# Save as TSV for easy inspection
write.table(univar_df, 
            file = file.path(output_dir, "PD_females_uacr_female.univ.ukb.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(bivar_df,
            file = file.path(output_dir, "PD_females_uacr_female.bivar.ukb.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# Print summary
print(paste("\nAnalysis complete! Results saved to", output_dir))
print(paste("Number of univariate results:", nrow(univar_df)))
print(paste("Number of bivariate results:", nrow(bivar_df)))

if (nrow(bivar_df) > 0) {
  print(paste("Significant bivariate results (p < 0.05):", sum(bivar_df$p < 0.05, na.rm = TRUE)))
}