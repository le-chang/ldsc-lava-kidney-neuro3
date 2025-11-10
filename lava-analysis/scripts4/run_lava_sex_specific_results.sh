#!/bin/bash
#SBATCH --job-name=lava_sex_spec
#SBATCH --account=def-gsarah
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=lava_sex_specific_%j.out
#SBATCH --error=lava_sex_specific_%j.err

# Description: Run sex-specific pairwise LAVA analyses
# Females: eGFR-PD, hematuria-PD, uACR-PD
# Males: eGFR-PD, hematuria-PD, uACR-PD

module load StdEnv/2023
module load r/4.3.1

# Set R library path
export R_LIBS=~/.local/R/4.3.1/:$R_LIBS

# Define paths
PROJECT_DIR="/home/lchang24/projects/def-gsarah/lchang24/github/ldsc-lava-kidney-neuro3/lava-analysis"
SCRIPT_DIR="${PROJECT_DIR}/scripts4"

echo "========================================="
echo "Sex-Specific LAVA Analysis"
echo "========================================="
echo "Start time: $(date)"
echo "Memory allocated: 64GB"
echo ""

Rscript /home/lchang24/projects/def-gsarah/lchang24/github/ldsc-lava-kidney-neuro3/lava-analysis/scripts4/extract_sex_specific_lava_results.R

