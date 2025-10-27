#!/bin/bash
# Install missing R packages for LAVA pipeline
# Run this directly in your terminal (NOT with sbatch)

echo "=========================================="
echo "Installing Missing R Packages"
echo "=========================================="
echo ""
echo "Missing packages to install:"
echo "  - tidyverse"
echo "  - stringr"
echo "  - dplyr"
echo "  - biglasso"
echo "  - remotes"
echo ""
echo "Already installed:"
echo "  ✓ LAVA (0.1.5)"
echo "  ✓ data.table (1.17.8)"
echo "  ✓ Matrix (1.5.4.1)"
echo ""
echo "This will take approximately 20-40 minutes..."
echo "Starting in 3 seconds (Press Ctrl+C to cancel)..."
sleep 3

# Set R library path
export R_LIBS=/home/lchang24/.local/R/4.3.1/

# Load R module
module load r/4.3.1

echo ""
echo "R Library Path: $R_LIBS"
echo "Starting installation..."
echo ""

# Install packages
Rscript -e "
# Set library path
.libPaths('/home/lchang24/.local/R/4.3.1/')

# Verify library path
cat('Library paths:\n')
print(.libPaths())
cat('\n')

# Set CRAN mirror
options(repos = c(CRAN = 'https://cloud.r-project.org'))

# Packages to install
missing_packages <- c('remotes', 'biglasso', 'tidyverse', 'stringr', 'dplyr')

# Install function
install_pkg <- function(pkg) {
  cat(sprintf('\\n==========================================\\n'))
  cat(sprintf('Installing: %s\\n', pkg))
  cat(sprintf('==========================================\\n'))
  
  tryCatch({
    install.packages(pkg, lib = '/home/lchang24/.local/R/4.3.1/', dependencies = TRUE)
    
    if (require(pkg, character.only = TRUE, quietly = TRUE)) {
      cat(sprintf('✓ SUCCESS: %s (version %s)\\n', pkg, packageVersion(pkg)))
      return(TRUE)
    } else {
      cat(sprintf('✗ FAILED: %s could not be loaded\\n', pkg))
      return(FALSE)
    }
  }, error = function(e) {
    cat(sprintf('✗ ERROR installing %s: %s\\n', pkg, e\$message))
    return(FALSE)
  })
}

# Install packages in order
# 1. Install remotes first (needed for GitHub packages)
cat('\\n##############################################\\n')
cat('STEP 1/5: Installing remotes\\n')
cat('##############################################\\n')
install_pkg('remotes')

# 2. Install biglasso
cat('\\n##############################################\\n')
cat('STEP 2/5: Installing biglasso\\n')
cat('##############################################\\n')
install_pkg('biglasso')

# 3. Install tidyverse (this takes the longest)
cat('\\n##############################################\\n')
cat('STEP 3/5: Installing tidyverse (this may take 15-30 min)\\n')
cat('##############################################\\n')
install_pkg('tidyverse')

# 4. Install stringr (might be included in tidyverse, but install separately to be sure)
cat('\\n##############################################\\n')
cat('STEP 4/5: Installing stringr\\n')
cat('##############################################\\n')
install_pkg('stringr')

# 5. Install dplyr (might be included in tidyverse, but install separately to be sure)
cat('\\n##############################################\\n')
cat('STEP 5/5: Installing dplyr\\n')
cat('##############################################\\n')
install_pkg('dplyr')

# Final verification
cat('\\n##############################################\\n')
cat('FINAL VERIFICATION\\n')
cat('##############################################\\n\\n')

required <- c('LAVA', 'tidyverse', 'stringr', 'dplyr', 'data.table', 
              'Matrix', 'biglasso', 'remotes')

all_ok <- TRUE
for (pkg in required) {
  if (require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(sprintf('✓ %-15s version %s\\n', pkg, packageVersion(pkg)))
  } else {
    cat(sprintf('✗ %-15s MISSING\\n', pkg))
    all_ok <- FALSE
  }
}

cat('\\n')
if (all_ok) {
  cat('========================================\\n')
  cat('✓ ALL PACKAGES INSTALLED SUCCESSFULLY!\\n')
  cat('========================================\\n')
  cat('\\nYou can now run your LAVA pipeline:\\n')
  cat('  sbatch lava_full_pipeline.sh\\n\\n')
} else {
  cat('========================================\\n')
  cat('✗ SOME PACKAGES STILL MISSING\\n')
  cat('========================================\\n')
  cat('\\nPlease review the errors above.\\n')
  cat('You may need to install failed packages manually in R.\\n\\n')
}
"

echo ""
echo "=========================================="
echo "Installation script completed!"
echo "=========================================="
echo ""
echo "Run 'bash check_all_r_pkgs.sh' to verify all packages are installed."
