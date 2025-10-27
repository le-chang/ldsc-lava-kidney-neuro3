#!/bin/bash
# Description: Install required R packages for LAVA analysis
# Run this directly on a login node (NOT through SLURM)
# Usage: bash install_packages.sh

echo "====================================="
echo "R Package Installer for LAVA Analysis"
echo "====================================="
echo "Date: $(date)"
echo ""

# Load R module
echo "Loading R 4.3.1..."
module load r/4.3.1

echo "R version:"
R --version | head -n 1
echo ""

# Check if R is loaded
if ! command -v R &> /dev/null; then
    echo "Error: R is not available. Please check module load."
    exit 1
fi

echo "====================================="
echo "Installing R packages..."
echo "====================================="
echo ""

# Run R script to install packages
R --vanilla --no-save << 'EOF'

cat("Starting package installation...\n\n")

# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

# List of required packages
required_packages <- c(
  "tidyverse",      # Data manipulation
  "data.table",     # Fast data reading/writing
  "stringr",        # String manipulation
  "gtools",         # General tools
  "here",           # Path management
  "readr",          # Reading rectangular data
  "devtools"        # Needed for LAVA installation
)

# Function to check if a package is installed
is_installed <- function(pkg) {
  return(pkg %in% installed.packages()[, "Package"])
}

# Install CRAN packages
cat("===================================\n")
cat("Installing CRAN packages...\n")
cat("===================================\n\n")

for (pkg in required_packages) {
  if (is_installed(pkg)) {
    cat(paste("✓", pkg, "- already installed\n"))
  } else {
    cat(paste("\n>>> Installing", pkg, "...\n"))
    tryCatch({
      install.packages(pkg, dependencies = TRUE)
      cat(paste("✓", pkg, "- SUCCESS\n"))
    }, error = function(e) {
      cat(paste("✗ Error installing", pkg, ":", e$message, "\n"))
    })
  }
}

# Install LAVA from GitHub
cat("\n===================================\n")
cat("Installing LAVA from GitHub...\n")
cat("===================================\n\n")

if (is_installed("LAVA")) {
  cat("✓ LAVA - already installed\n")
} else {
  cat(">>> Installing LAVA from GitHub (josefin-werme/LAVA)...\n")
  tryCatch({
    library(devtools)
    devtools::install_github("josefin-werme/LAVA", upgrade = "never")
    cat("✓ LAVA - SUCCESS\n")
  }, error = function(e) {
    cat(paste("✗ Error installing LAVA:", e$message, "\n"))
    cat("\nIf installation failed, you may need to:\n")
    cat("1. Check internet connectivity\n")
    cat("2. Ensure you have a C++ compiler available\n")
    cat("3. Try manual installation in R console\n")
  })
}

# Final verification
cat("\n===================================\n")
cat("VERIFICATION\n")
cat("===================================\n\n")

all_packages <- c(required_packages, "LAVA")
installed_count <- 0

for (pkg in all_packages) {
  if (is_installed(pkg)) {
    version <- tryCatch(
      as.character(packageVersion(pkg)),
      error = function(e) "unknown"
    )
    cat(paste("✓", pkg, "version", version, "\n"))
    installed_count <- installed_count + 1
  } else {
    cat(paste("✗", pkg, "- NOT INSTALLED\n"))
  }
}

cat("\n===================================\n")
cat("SUMMARY\n")
cat("===================================\n")
cat(paste("Successfully installed:", installed_count, "of", length(all_packages), "packages\n\n"))

if (installed_count == length(all_packages)) {
  cat("✓ SUCCESS! All packages are ready.\n")
  cat("You can now run the LAVA analysis.\n")
} else {
  cat("✗ WARNING! Some packages failed to install.\n")
  cat("Please install them manually.\n")
}

EOF

echo ""
echo "====================================="
echo "Installation script completed"
echo "====================================="
echo ""
echo "Next steps:"
echo "1. Check the output above for any errors"
echo "2. If all packages installed successfully, run:"
echo "   sbatch run_lava_comparison.sh"
echo ""
