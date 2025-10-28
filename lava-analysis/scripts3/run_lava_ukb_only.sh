#!/bin/bash
#SBATCH --job-name=lava_ukb
#SBATCH --account=def-gsarah
#SBATCH --time=20:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --output=lava_ukb_%j.out
#SBATCH --error=lava_ukb_%j.err

# Description: Run LAVA analysis with UKB reference only
# This allows running references separately if needed

module load StdEnv/2023
module load r/4.3.1

# Set R library path
export R_LIBS=~/.local/R/4.3.1/:$R_LIBS

# Define paths
PROJECT_DIR="/home/lchang24/projects/def-gsarah/lchang24/github/ldsc-lava-kidney-neuro3/lava-analysis"
SCRIPT_DIR="${PROJECT_DIR}/scripts3"

echo "========================================="
echo "LAVA Analysis - UKB Reference Only"
echo "========================================="
echo "Start time: $(date)"
echo ""

# Run analysis
Rscript ${SCRIPT_DIR}/lava_analysis_ukb.R

echo ""
echo "========================================="
echo "UKB Analysis Complete"
echo "========================================="
echo "End time: $(date)"
