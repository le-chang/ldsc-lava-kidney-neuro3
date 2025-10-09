#!/bin/bash


#SBATCH --time=20:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --job-name=diagnose
#SBATCH --output=logs/diagnose%A_%a.out
#SBATCH --error=logs/diagnose%A_%a.err
#SBATCH --account=def-gsarah

set -x  # Print each command

module load StdEnv/2023
module load apptainer/1.3.5

SIF_FILE="/home/lchang24/projects/def-gsarah/lchang24/github/nf-genetic-correlations_fork/nf-genetic-correlations/bin/ldsc_latest.sif"

# Test 1: Can we run apptainer?
echo "Test 1: Basic apptainer test"
apptainer exec ${SIF_FILE} python --version

# Test 2: Can munge_sumstats be found?
echo "Test 2: Testing munge_sumstats"
apptainer exec ${SIF_FILE} which munge_sumstats.py

# Test 3: Run munge_sumstats with --help
echo "Test 3: munge_sumstats help"
apptainer exec ${SIF_FILE} munge_sumstats.py --help

# Test 4: Try ONE munge with minimal data
echo "Test 4: Running one munge command"
apptainer exec ${SIF_FILE} munge_sumstats.py \
    --sumstats /home/lchang24/projects/def-gsarah/lchang24/github/nf-genetic-correlations_fork/nf-genetic-correlations/data/sumstats_original/AD_wightman2021_with_rsid.txt \
    --out /home/lchang24/scratch/test_output \
    --N 762917 \
    --snp SNP \
    --a1 A1 \
    --a2 A2 \
    --p P \
    --signed-sumstats BETA,0

echo "Diagnostic complete!"
