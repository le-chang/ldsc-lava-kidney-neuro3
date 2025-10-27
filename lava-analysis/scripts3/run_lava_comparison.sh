#!/bin/bash
#SBATCH --job-name=lava_1000g_ukb_comparison
#SBATCH --account=def-gsarah
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=lava_comparison_%j.out
#SBATCH --error=lava_comparison_%j.err

# Description: Run LAVA analysis comparing 1000G vs UKB reference
# for PD_females vs uacr_female trait pair

module load StdEnv/2023
module load r/4.3.1

export R_LIBS=~/.local/R/4.3.1/

# Define paths
PROJECT_DIR="/home/lchang24/projects/def-gsarah/lchang24/github/ldsc-lava-kidney-neuro3/lava-analysis"
SCRIPT_DIR="${PROJECT_DIR}/scripts3"
OUTPUT_DIR="${PROJECT_DIR}/results"

# Create output directory if it doesn't exist
mkdir -p ${OUTPUT_DIR}

echo "Starting LAVA analysis comparison..."
echo "Using R 4.3.1"
echo "Project directory: ${PROJECT_DIR}"
echo ""

# Run LAVA analysis with 1000G reference
echo "Running LAVA with 1000G reference..."
Rscript ${SCRIPT_DIR}/lava_analysis_1000g.R

# Run LAVA analysis with UKB reference
echo "Running LAVA with UKB reference..."
Rscript ${SCRIPT_DIR}/lava_analysis_ukb.R

# Compare results
echo "Comparing results..."
Rscript ${SCRIPT_DIR}/compare_lava_results.R

echo "Analysis complete!"