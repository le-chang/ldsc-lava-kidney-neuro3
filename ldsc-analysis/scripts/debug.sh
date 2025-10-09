#!/bin/bash
#SBATCH --time=0:30:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --job-name=debug_lava
#SBATCH --output=logs/debug_lava_%j.out
#SBATCH --error=logs/debug_lava_%j.err
#SBATCH --account=def-gsarah

set -x  # Print every command

echo "Starting debug at: $(date)"

module load StdEnv/2023
module load apptainer/1.3.5

SIF_FILE="/home/lchang24/projects/def-gsarah/lchang24/github/nf-genetic-correlations_fork/nf-genetic-correlations/bin/ldsc_latest.sif"

# Test 1: Can we read the .lava.gz file?
echo "Test 1: Checking file"
ls -lh /home/lchang24/projects/def-gsarah/sumstats/kidney_neurodegen/PDnp_females.lava.gz

# Test 2: Can we decompress it?
echo "Test 2: Testing decompression"
zcat /home/lchang24/projects/def-gsarah/sumstats/kidney_neurodegen/PDnp_females.lava.gz | head -5

# Test 3: Check what columns it has
echo "Test 3: File header"
zcat /home/lchang24/projects/def-gsarah/sumstats/kidney_neurodegen/PDnp_females.lava.gz | head -1

# Test 4: Test apptainer with explicit binds
echo "Test 4: Testing apptainer with binds"
apptainer exec \
    --bind /home/lchang24:/home/lchang24 \
    --bind /scratch:/scratch \
    ${SIF_FILE} \
    python --version

# Test 5: Try the actual command with verbose output
echo "Test 5: Running munge with smallest file"
apptainer exec \
    --bind /home/lchang24:/home/lchang24 \
    --bind /scratch:/scratch \
    ${SIF_FILE} \
    munge_sumstats.py \
    --sumstats /home/lchang24/projects/def-gsarah/sumstats/kidney_neurodegen/PDnp_females.lava.gz \
    --out /home/lchang24/scratch/test_pdnp_females \
    --N 19773 \
    --ignore N \
    --snp SNP \
    --a1 A1 \
    --a2 A2 \
    --p P \
    --signed-sumstats BETA,0

echo "Debug complete at: $(date)"
