#!/bin/bash
# Script to check and install LAVA package

# Set R library path
export R_LIBS=~/.local/R/4.3.1/

# Load R module
module load r/4.3.1

echo "=========================================="
echo "Checking LAVA Installation"
echo "=========================================="
echo ""
echo "R Library Path: $R_LIBS"
echo ""

# Check if LAVA is installed
Rscript -e "
if (require('LAVA', quietly = TRUE)) {
  cat('✓ LAVA is installed!\n')
  cat('Version:', as.character(packageVersion('LAVA')), '\n')
  cat('Location:', find.package('LAVA'), '\n')
} else {
  cat('✗ LAVA is NOT installed.\n')
  cat('\nWould you like to install it? (This script can do it automatically)\n')
}
"

echo ""
echo "=========================================="
echo "To install LAVA if not present, run:"
echo "sbatch install_lava.sh"
echo "=========================================="
