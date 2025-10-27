# Description: Prepare input.info.txt file for LAVA analysis
# This creates the input info file needed by LAVA from your summary statistics

library(tidyverse)

# Set paths
project_dir <- "/home/lchang24/projects/def-gsarah"
sumstats_dir <- file.path(project_dir, "sumstats/kidney_neurodegen")

# Define your summary statistics files
sumstats_files <- data.frame(
  phenotype = c("PD_females", "uacr_female"),
  filename = c(
    file.path(sumstats_dir, "PD_females.lava.gz"),
    "/home/lchang24/projects/def-gsarah/lchang24/github/nf-genetic-correlations_fork/nf-genetic-correlations/data/sumstats_original/uacr_female.tsv"
  ),
  # Sample sizes from your data
  cases = c(NA, NA),  # Set to NA for continuous traits, or actual case count for binary
  controls = c(NA, NA),  # Set to NA for continuous traits, or actual control count for binary
  prevalence = c(NA, NA)  # Set to NA for continuous traits, or prevalence for binary
)

# Function to detect column names and create input info
prepare_lava_input <- function(pheno_name, file_path) {
  
  cat(sprintf("Processing %s...\n", pheno_name))
  
  # Read header to determine column names
  if (grepl("\\.gz$", file_path)) {
    header <- readLines(gzfile(file_path), n = 1)
  } else {
    header <- readLines(file_path, n = 1)
  }
  
  cols <- strsplit(header, "\t")[[1]]
  
  # Map column names (adjust these mappings based on your actual column names)
  col_mapping <- list(
    SNP = ifelse("SNP" %in% cols, "SNP", 
                 ifelse("variant_ID" %in% cols, "variant_ID", NA)),
    A1 = ifelse("A1" %in% cols, "A1", NA),
    A2 = ifelse("A2" %in% cols, "A2", NA),
    N = ifelse("N" %in% cols, "N", NA),
    Z = NA,  # Will calculate from BETA/SE or use P-value
    BETA = ifelse("BETA" %in% cols, "BETA", NA),
    SE = ifelse("SE" %in% cols, "SE", NA),
    P = ifelse("P" %in% cols, "P", NA)
  )
  
  return(data.frame(
    phenotype = pheno_name,
    filename = file_path,
    SNP = col_mapping$SNP,
    A1 = col_mapping$A1,
    A2 = col_mapping$A2,
    N = col_mapping$N,
    BETA = col_mapping$BETA,
    SE = col_mapping$SE,
    P = col_mapping$P
  ))
}

# Create input info for each phenotype
input_info <- lapply(1:nrow(sumstats_files), function(i) {
  prepare_lava_input(
    sumstats_files$phenotype[i],
    sumstats_files$filename[i]
  )
}) %>% bind_rows()

# Write input.info.txt file
output_file <- file.path(project_dir, "info_files", "input.info.txt")
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

write.table(
  input_info,
  file = output_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat(sprintf("\nInput info file created: %s\n", output_file))
cat("\nPlease review the file and adjust column mappings if needed!\n")
cat("Also ensure you have created the sample_overlap.txt file.\n")

# Print what was created
print(input_info)

# Example sample_overlap.txt structure (you need to create this separately)
cat("\n=== Sample overlap file example ===\n")
cat("You need to create a sample_overlap.txt file with the following structure:\n")
cat("phenotype1\tphenotype2\toverlap\n")
cat("PD_females\tuacr_female\t0\n")
cat("\nWhere overlap is the proportion of sample overlap between the two traits (0 to 1)\n")
