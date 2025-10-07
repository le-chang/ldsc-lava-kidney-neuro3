#!/bin/bash

#SBATCH --time=4:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --job-name=munge_lava_fix
#SBATCH --output=logs/munge_lava_%A_%a.out
#SBATCH --error=logs/munge_lava_%A_%a.err
#SBATCH --account=def-gsarah

module load StdEnv/2023
module load apptainer/1.3.5

SIF_FILE="/home/lchang24/projects/def-gsarah/lchang24/github/nf-genetic-correlations_fork/nf-genetic-correlations/bin/ldsc_latest.sif"
OUTPUT_DIR="/home/lchang24/scratch/ldsc_munged"
mkdir -p ${OUTPUT_DIR}

# ============================================================================
# FIX: Add --ignore N to tell LDSC to ignore the N column from the file
# ============================================================================

# ============================================================================
# PD_females.lava.gz - FIXED
# ============================================================================
echo "Processing PD_females - normal orientation..."
apptainer run -W ${SLURM_TMPDIR} ${SIF_FILE} munge_sumstats.py \
    --sumstats /home/lchang24/projects/def-gsarah/sumstats/kidney_neurodegen/PD_females.lava.gz \
    --out ${OUTPUT_DIR}/PD_females_normal \
    --N 104082 \
    --ignore N \
    --snp SNP \
    --a1 A1 \
    --a2 A2 \
    --p P \
    --signed-sumstats BETA,0

echo "Processing PD_females - flipped orientation..."
apptainer run -W ${SLURM_TMPDIR} ${SIF_FILE} munge_sumstats.py \
    --sumstats /home/lchang24/projects/def-gsarah/sumstats/kidney_neurodegen/PD_females.lava.gz \
    --out ${OUTPUT_DIR}/PD_females_flipped \
    --N 104082 \
    --ignore N \
    --snp SNP \
    --a1 A2 \
    --a2 A1 \
    --p P \
    --signed-sumstats BETA,0

# ============================================================================
# PD_males.lava.gz - FIXED
# ============================================================================
echo "Processing PD_males - normal orientation..."
apptainer run -W ${SLURM_TMPDIR} ${SIF_FILE} munge_sumstats.py \
    --sumstats /home/lchang24/projects/def-gsarah/sumstats/kidney_neurodegen/PD_males.lava.gz \
    --out ${OUTPUT_DIR}/PD_males_normal \
    --N 110616 \
    --ignore N \
    --snp SNP \
    --a1 A1 \
    --a2 A2 \
    --p P \
    --signed-sumstats BETA,0

echo "Processing PD_males - flipped orientation..."
apptainer run -W ${SLURM_TMPDIR} ${SIF_FILE} munge_sumstats.py \
    --sumstats /home/lchang24/projects/def-gsarah/sumstats/kidney_neurodegen/PD_males.lava.gz \
    --out ${OUTPUT_DIR}/PD_males_flipped \
    --N 110616 \
    --ignore N \
    --snp SNP \
    --a1 A2 \
    --a2 A1 \
    --p P \
    --signed-sumstats BETA,0

# ============================================================================
# PDnp_females.lava.gz - FIXED
# ============================================================================
echo "Processing PDnp_females - normal orientation..."
apptainer run -W ${SLURM_TMPDIR} ${SIF_FILE} munge_sumstats.py \
    --sumstats /home/lchang24/projects/def-gsarah/sumstats/kidney_neurodegen/PDnp_females.lava.gz \
    --out ${OUTPUT_DIR}/PDnp_females_normal \
    --N 19773 \
    --ignore N \
    --snp SNP \
    --a1 A1 \
    --a2 A2 \
    --p P \
    --signed-sumstats BETA,0

echo "Processing PDnp_females - flipped orientation..."
apptainer run -W ${SLURM_TMPDIR} ${SIF_FILE} munge_sumstats.py \
    --sumstats /home/lchang24/projects/def-gsarah/sumstats/kidney_neurodegen/PDnp_females.lava.gz \
    --out ${OUTPUT_DIR}/PDnp_females_flipped \
    --N 19773 \
    --ignore N \
    --snp SNP \
    --a1 A2 \
    --a2 A1 \
    --p P \
    --signed-sumstats BETA,0

# ============================================================================
# PDnp_males.lava.gz - FIXED
# ============================================================================
echo "Processing PDnp_males - normal orientation..."
apptainer run -W ${SLURM_TMPDIR} ${SIF_FILE} munge_sumstats.py \
    --sumstats /home/lchang24/projects/def-gsarah/sumstats/kidney_neurodegen/PDnp_males.lava.gz \
    --out ${OUTPUT_DIR}/PDnp_males_normal \
    --N 24053 \
    --ignore N \
    --snp SNP \
    --a1 A1 \
    --a2 A2 \
    --p P \
    --signed-sumstats BETA,0

echo "Processing PDnp_males - flipped orientation..."
apptainer run -W ${SLURM_TMPDIR} ${SIF_FILE} munge_sumstats.py \
    --sumstats /home/lchang24/projects/def-gsarah/sumstats/kidney_neurodegen/PDnp_males.lava.gz \
    --out ${OUTPUT_DIR}/PDnp_males_flipped \
    --N 24053 \
    --ignore N \
    --snp SNP \
    --a1 A2 \
    --a2 A1 \
    --p P \
    --signed-sumstats BETA,0

echo ""
echo "=========================================="
echo "All .lava.gz files processed successfully!"
echo "End Time: $(date)"
echo "=========================================="
