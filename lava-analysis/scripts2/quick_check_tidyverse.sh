#!/bin/bash
# Quick one-line check for tidyverse

export R_LIBS=~/.local/R/4.3.1/
module load r/4.3.1

echo "Checking tidyverse..."
Rscript -e ".libPaths('~/.local/R/4.3.1/'); if (require('tidyverse', quietly=TRUE)) { cat('✓ tidyverse is installed (version', as.character(packageVersion('tidyverse')), ')\n') } else { cat('✗ tidyverse is NOT installed\n') }"
