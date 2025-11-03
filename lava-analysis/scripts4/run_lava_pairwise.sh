#!/bin/bash
#SBATCH --job-name=lava_pairwise
#SBATCH --account=def-gsarah
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=lava_pairwise_%j.out
#SBATCH --error=lava_pairwise_%j.err

# Description: Run pairwise LAVA analyses
# 1. hematuria vs PD
# 2. eGFR vs PD
# 3. uacr vs PD

module load StdEnv/2023
module load r/4.3.1

# Set R library path
export R_LIBS=~/.local/R/4.3.1/:$R_LIBS

# Define paths
PROJECT_DIR="/home/lchang24/projects/def-gsarah/lchang24/github/ldsc-lava-kidney-neuro3/lava-analysis"
SCRIPT_DIR="${PROJECT_DIR}/scripts4"

echo "========================================="
echo "LAVA Pairwise Analyses"
echo "========================================="
echo "Start time: $(date)"
echo ""
echo "Running 3 pairwise analyses:"
echo "1. hematuria vs PD"
echo "2. eGFR vs PD"
echo "3. uacr vs PD"
echo ""

# Run pairwise analyses
Rscript ${SCRIPT_DIR}/lava_pairwise_analyses.R

if [ $? -eq 0 ]; then
    echo ""
    echo "========================================="
    echo "✓ All pairwise analyses completed successfully"
    echo "========================================="
else
    echo ""
    echo "========================================="
    echo "✗ Error in pairwise analyses"
    echo "========================================="
    exit 1
fi

echo ""
echo "End time: $(date)"
