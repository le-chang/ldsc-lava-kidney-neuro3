#!/bin/bash
#SBATCH --job-name=lava_4traits
#SBATCH --account=def-gsarah
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=lava_4traits_%j.out
#SBATCH --error=lava_4traits_%j.err

# Description: Run LAVA analysis for 4 traits
# uacr, eGFR, hematuria, and PD

module load StdEnv/2023
module load r/4.3.1

# Set R library path
export R_LIBS=~/.local/R/4.3.1/:$R_LIBS

# Define paths
PROJECT_DIR="/home/lchang24/projects/def-gsarah/lchang24/github/ldsc-lava-kidney-neuro3/lava-analysis"
SCRIPT_DIR="${PROJECT_DIR}/scripts4"

echo "========================================="
echo "LAVA Analysis - 4 Traits"
echo "========================================="
echo "Start time: $(date)"
echo ""

# Step 1: Prepare summary statistics
echo "Step 1: Preparing summary statistics..."
Rscript ${SCRIPT_DIR}/prep_sumstats_4traits.R

if [ $? -eq 0 ]; then
    echo "✓ Summary statistics preparation completed successfully"
else
    echo "✗ Error in summary statistics preparation"
    exit 1
fi

echo ""

# Step 2: Run LAVA analysis
echo "Step 2: Running LAVA analysis..."
Rscript ${SCRIPT_DIR}/lava_analysis_4traits.R

if [ $? -eq 0 ]; then
    echo "✓ LAVA analysis completed successfully"
else
    echo "✗ Error in LAVA analysis"
    exit 1
fi

echo ""
echo "========================================="
echo "4-Trait Analysis Complete"
echo "========================================="
echo "End time: $(date)"
