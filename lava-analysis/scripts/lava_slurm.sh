#!/bin/bash
#SBATCH --job-name=lava_ref_comparison
#SBATCH --account=def-gsarah
#SBATCH --time=12:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --output=lava_comparison_%j.out
#SBATCH --error=lava_comparison_%j.err

# Load required modules
module load r/4.3.1

# Set working directory
cd /home/lchang24/projects/def-gsarah/sumstats/kidney_neurodegen

# Create output directories if they don't exist
mkdir -p lava_results/1000G_ref
mkdir -p lava_results/UKB_ref

# Run LAVA with 1000G reference
echo "Running LAVA with 1000G reference..."
Rscript lava_comparison_1000G.R

# Run LAVA with UKB reference  
echo "Running LAVA with UKB reference..."
Rscript lava_comparison_UKB.R

echo "LAVA comparison completed!"
