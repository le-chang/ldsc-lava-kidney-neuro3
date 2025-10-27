#!/bin/bash

#SBATCH --job-name=ldsc_rg_ad_uacr
#SBATCH --account=def-gsarah
#SBATCH --time=3:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --array=0-2
#SBATCH --output=ldsc_rg_ad_uacr-%A-%a.out
#SBATCH --error=ldsc_rg_ad_uacr-%A-%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=le.chang@umontreal.ca

# Exit on error and print commands
set -e
set -x

echo "=========================================="
echo "LDSC rg Job started at: $(date)"
echo "Job ID: ${SLURM_ARRAY_JOB_ID}"
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Node: ${SLURM_NODELIST}"
echo "=========================================="
echo ""

# ============================================
# CONFIGURATION
# ============================================

module load StdEnv/2023
module load apptainer/1.3.5

SIF_FILE="/home/lchang24/projects/def-gsarah/lchang24/github/nf-genetic-correlations_fork/nf-genetic-correlations/bin/ldsc_latest.sif"
LD_SCORES_DIR="/home/lchang24/projects/def-gsarah/lchang24/github/nf-genetic-correlations_fork/nf-genetic-correlations/data/ld_reference/eur_w_ld_chr"
AD_FILE="/home/lchang24/projects/def-gsarah/sumstats/AD_wightman2021.sumstats.gz"
UACR_DIR="/home/lchang24/scratch/munge_output_additional"
OUTPUT_DIR="/home/lchang24/scratch/ldsc_rg_results_ad_uacr"

# Create output directory if it doesn't exist
mkdir -p ${OUTPUT_DIR}

# Define UACR traits (3 phenotypes, using A1ref versions)
UACR_FILES=(
    "${UACR_DIR}/uacr_male.tsv.A1ref.munged.sumstats.gz"
    "${UACR_DIR}/uacr_sex_combined.tsv.A1ref.munged.sumstats.gz"
    "${UACR_DIR}/uacr_female.tsv.A1ref.munged.sumstats.gz"
)

UACR_NAMES=(
    "uacr_male"
    "uacr_sex_combined"
    "uacr_female"
)

# Get the specific UACR trait for this task
UACR_FILE="${UACR_FILES[$SLURM_ARRAY_TASK_ID]}"
UACR_NAME="${UACR_NAMES[$SLURM_ARRAY_TASK_ID]}"

echo "=========================================="
echo "Processing genetic correlation:"
echo "AD trait: AD_wightman2021"
echo "UACR trait: ${UACR_NAME}"
echo "=========================================="
echo ""

# Check if input files exist
if [ ! -f "${AD_FILE}" ]; then
    echo "ERROR: AD file not found: ${AD_FILE}"
    exit 1
fi

if [ ! -f "${UACR_FILE}" ]; then
    echo "ERROR: UACR file not found: ${UACR_FILE}"
    exit 1
fi

# ============================================
# Run LDSC genetic correlation
# ============================================

OUTPUT_PREFIX="${OUTPUT_DIR}/rg_AD_wightman2021_vs_${UACR_NAME}"

echo "Running LDSC genetic correlation..."
echo "Started at: $(date)"
echo ""

time apptainer exec \
    --bind /home/lchang24:/home/lchang24 \
    --bind /scratch:/scratch \
    ${SIF_FILE} \
    ldsc.py \
    --rg ${AD_FILE},${UACR_FILE} \
    --ref-ld-chr ${LD_SCORES_DIR}/ \
    --w-ld-chr ${LD_SCORES_DIR}/ \
    --out ${OUTPUT_PREFIX}

echo ""
echo "================================================"
echo "COMPLETED!"
echo "================================================"
echo "Output: ${OUTPUT_PREFIX}.log"
echo "Finished at: $(date)"
echo ""

# Display results
if [ -f "${OUTPUT_PREFIX}.log" ]; then
    echo "===== LDSC rg Results ====="
    echo ""
    # Extract genetic correlation results
    grep -A 10 "Genetic Correlation" ${OUTPUT_PREFIX}.log || echo "Results processing..."
    echo ""
    echo "Full log file: ${OUTPUT_PREFIX}.log"
else
    echo "Warning: Output log file not found!"
    exit 1
fi

echo ""
echo "=========================================="
echo "Array task ${SLURM_ARRAY_TASK_ID} completed at: $(date)"
echo "AD: AD_wightman2021"
echo "UACR: ${UACR_NAME}"
echo "=========================================="
