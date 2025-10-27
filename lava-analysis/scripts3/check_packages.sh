#!/bin/bash
# Description: Check which R packages are installed (no installation)
# Run this on a login node to check package status
# Usage: bash check_packages.sh

echo "====================================="
echo "R Package Status Checker"
echo "====================================="
echo ""

# Load R module
module load r/4.3.1

echo "R version:"
R --version | head -n 1
echo ""

# Check packages
R --vanilla --no-save << 'EOF'

# List of required packages
required_packages <- c(
  "LAVA",
  "tidyverse",
  "data.table",
  "stringr",
  "gtools",
  "here",
  "readr",
  "devtools",
  "dplyr",
  "ggplot2"
)

cat("Checking installed packages...\n")
cat("================================\n\n")

installed <- c()
missing <- c()

for (pkg in required_packages) {
  if (pkg %in% installed.packages()[, "Package"]) {
    version <- as.character(packageVersion(pkg))
    cat(paste("✓", pkg, "version", version, "\n"))
    installed <- c(installed, pkg)
  } else {
    cat(paste("✗", pkg, "- NOT INSTALLED\n"))
    missing <- c(missing, pkg)
  }
}

cat("\n================================\n")
cat("SUMMARY\n")
cat("================================\n")
cat(paste("Installed:", length(installed), "\n"))
cat(paste("Missing:", length(missing), "\n"))

if (length(missing) > 0) {
  cat("\nMissing packages:\n")
  for (pkg in missing) {
    cat(paste("  -", pkg, "\n"))
  }
  cat("\nTo install missing packages, run:\n")
  cat("  bash install_packages.sh\n")
} else {
  cat("\n✓ All required packages are installed!\n")
}

EOF
