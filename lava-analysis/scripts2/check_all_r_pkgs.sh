#!/bin/bash
# Script to check all required R packages for LAVA pipeline

# Set R library path
export R_LIBS=~/.local/R/4.3.1/

# Load R module
module load r/4.3.1

echo "=========================================="
echo "Checking R Package Installation"
echo "=========================================="
echo ""
echo "R Library Path: $R_LIBS"
echo ""

# Check all required packages
Rscript -e "
# Set library path
.libPaths('~/.local/R/4.3.1/')

# List of required packages
required_packages <- c(
  'LAVA',
  'tidyverse',
  'stringr',
  'dplyr',
  'data.table',
  'Matrix',
  'biglasso',
  'remotes'
)

cat('Checking required packages:\n')
cat('==========================================\n\n')

# Check each package
all_installed <- TRUE
missing_packages <- c()

for (pkg in required_packages) {
  if (require(pkg, character.only = TRUE, quietly = TRUE)) {
    version <- as.character(packageVersion(pkg))
    cat(sprintf('✓ %-15s installed (version %s)\n', pkg, version))
  } else {
    cat(sprintf('✗ %-15s NOT INSTALLED\n', pkg))
    all_installed <- FALSE
    missing_packages <- c(missing_packages, pkg)
  }
}

cat('\n==========================================\n')

if (all_installed) {
  cat('✓ All required packages are installed!\n')
  cat('You can proceed with the LAVA analysis.\n')
} else {
  cat('✗ Missing packages:\n')
  for (pkg in missing_packages) {
    cat(sprintf('  - %s\n', pkg))
  }
  cat('\nTo install missing packages, run:\n')
  cat('  sbatch install_all_packages.sh\n')
}
"

echo ""
echo "=========================================="
