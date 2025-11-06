#!/bin/bash
#SBATCH --job-name=lava_sex_spec
#SBATCH --account=def-gsarah
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
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

# Step 1: Prepare summary statistics
echo "Step 1: Formatting sex-specific summary statistics..."
Rscript ${SCRIPT_DIR}/prep_sumstats_sex_specific.R

if [ $? -eq 0 ]; then
    echo "✓ Formatting completed successfully"
else
    echo "✗ Error in formatting"
    exit 1
fi

echo ""

# Step 2: Run LAVA analyses
echo "Step 2: Running sex-specific LAVA analyses..."
Rscript ${SCRIPT_DIR}/lava_pairwise_sex_specific.R

if [ $? -eq 0 ]; then
    echo "✓ LAVA analyses completed successfully"
else
    echo "✗ Error in LAVA analyses"
    exit 1
fi

echo ""
echo "========================================="
echo "Sex-Specific Analysis Complete"
echo "========================================="
echo "End time: $(date)"