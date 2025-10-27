#!/bin/bash
#SBATCH --job-name=lava_full_pipeline
#SBATCH --account=def-gsarah
#SBATCH --time=12:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --output=lava_pipeline_%j.out
#SBATCH --error=lava_pipeline_%j.err

module load StdEnv/2023
module load java/21.0.1
module load nextflow/24.10.2
module load r/4.3.1
module load apptainer/1.3.5

echo "=========================================="
echo "LAVA Analysis Pipeline"
echo "Starting at: $(date)"
echo "=========================================="

# Set R library path
export R_LIBS=~/.local/R/4.3.1/
echo "R Library Path: $R_LIBS"

# Load required modules
module load r/4.3.1

# Set working directory
cd /home/lchang24/projects/def-gsarah/lchang24/github/ldsc-lava-kidney-neuro3/lava-analysis/scripts2
# cd /home/lchang24/projects/def-gsarah/sumstats/kidney_neurodegen

# Create output directories if they don't exist
mkdir -p lava_results/1000G_ref
mkdir -p lava_results/UKB_ref

# Step 1: Prepare LAVA input files
echo ""
echo "=========================================="
echo "Step 1: Preparing LAVA input files..."
echo "=========================================="
Rscript prepare_lava_input.R

if [ $? -ne 0 ]; then
    echo "ERROR: prepare_lava_input.R failed!"
    exit 1
fi

echo "Input preparation completed successfully!"

# Step 2: Run LAVA with 1000G reference
echo ""
echo "=========================================="
echo "Step 2: Running LAVA with 1000G reference..."
echo "=========================================="
Rscript lava_comparison_1000G.R

if [ $? -ne 0 ]; then
    echo "ERROR: lava_comparison_1000G.R failed!"
    exit 1
fi

echo "1000G analysis completed successfully!"

# Step 3: Run LAVA with UKB reference  
echo ""
echo "=========================================="
echo "Step 3: Running LAVA with UKB reference..."
echo "=========================================="
Rscript lava_comparison_UKB.R

if [ $? -ne 0 ]; then
    echo "ERROR: lava_comparison_UKB.R failed!"
    exit 1
fi

echo "UKB analysis completed successfully!"

# Summary
echo ""
echo "=========================================="
echo "Pipeline completed successfully!"
echo "Finished at: $(date)"
echo "=========================================="
echo ""
echo "Results saved to:"
echo "  - 1000G reference: lava_results/1000G_ref/"
echo "  - UKB reference: lava_results/UKB_ref/"
echo ""
echo "Output files:"
ls -lh lava_results/1000G_ref/
ls -lh lava_results/UKB_ref/
