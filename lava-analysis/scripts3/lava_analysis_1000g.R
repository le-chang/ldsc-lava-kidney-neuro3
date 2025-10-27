# Description: Run LAVA analysis with 1000G reference
# for PD_females vs uacr_female trait pair

# Load packages -----------------------------------------------------------
library(LAVA)
library(tidyverse)
library(data.table)

# Set arguments -----------------------------------------------------------
project_dir <- "/home/lchang24/projects/def-gsarah/lchang24/github/ldsc-lava-kidney-neuro3/lava-analysis"

# Define paths
ref_prefix <- "/home/lchang24/projects/def-gsarah/lchang24/github/nf-genetic-correlations_fork/nf-genetic-correlations/data/ld_reference/lava_ld_ref/g1000_eur/g1000_eur"
loc_file <- "/home/lchang24/projects/def-gsarah/lchang24/github/nf-genetic-correlations_fork/nf-genetic-correlations/data/lava_locus_file/blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile"
sumstats_dir <- "/home/lchang24/projects/def-gsarah/sumstats/kidney_neurodegen"
output_dir <- file.path(project_dir, "results", "1000g")

# Create output directory
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Prepare input info file -------------------------------------------------
# Create info file for LAVA
info_data <- data.frame(
  phenotype = c("PD_females", "uacr_female"),
  cases = c(NA, NA),  # Adjust if you have case counts
  controls = c(NA, NA),  # Adjust if you have control counts
  filename = c(
    file.path(sumstats_dir, "formatted_PD_females.tsv"),
    "/home/lchang24/projects/def-gsarah/lchang24/github/nf-genetic-correlations_fork/nf-genetic-correlations/data/sumstats_original/formatted_uacr_female.tsv"
  )
)

info_file <- file.path(project_dir, "scripts3", "info_file_1000g.txt")
write.table(info_data, info_file, sep = "\t", row.names = FALSE, quote = FALSE)

# Create sample overlap file ----------------------------------------------
# For now, creating a simple 2x2 identity matrix (assuming no overlap)
# Adjust based on your actual sample overlap information
sample_overlap <- matrix(c(1, 0, 0, 1), nrow = 2, ncol = 2)
rownames(sample_overlap) <- colnames(sample_overlap) <- c("PD_females", "uacr_female")

sample_overlap_file <- file.path(project_dir, "scripts3", "sample_overlap_1000g.txt")
write.table(sample_overlap, sample_overlap_file, sep = "\t", quote = FALSE)

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
print(paste("Starting LAVA analysis with 1000G reference for", n_loci, "loci"))

# Print progress
progress <- quantile(x = 1:n_loci, probs = seq(.05, 1, .05)) %>% ceiling()

# Set univariate threshold
univar_threshold <- 0.05 / n_loci

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
saveRDS(univar, file = file.path(output_dir, "PD_females_uacr_female.univ.1000g.rds"))
saveRDS(bivar, file = file.path(output_dir, "PD_females_uacr_female.bivar.1000g.rds"))

# Save as TSV for easy inspection
write.table(univar_df, 
            file = file.path(output_dir, "PD_females_uacr_female.univ.1000g.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(bivar_df,
            file = file.path(output_dir, "PD_females_uacr_female.bivar.1000g.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# Print summary
print(paste("Analysis complete! Results saved to", output_dir))
print(paste("Number of univariate results:", nrow(univar_df)))
print(paste("Number of bivariate results:", nrow(bivar_df)))