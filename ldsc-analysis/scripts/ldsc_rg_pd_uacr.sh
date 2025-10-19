#!/bin/bash

#SBATCH --job-name=ldsc_rg
#SBATCH --account=def-gsarah
#SBATCH --time=20:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --array=0-83
#SBATCH --output=ldsc_rg-%A-%a.out
#SBATCH --error=ldsc_rg-%A-%a.err
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
PD_DIR="/home/lchang24/scratch/munge_output"
UACR_DIR="/home/lchang24/scratch/munge_output_additional"
PD_ADDITIONAL_DIR="/home/lchang24/scratch/munge_output_additional"
OUTPUT_DIR="/home/lchang24/scratch/ldsc_rg_results"

# Create output directory if it doesn't exist
mkdir -p ${OUTPUT_DIR}

# Define ALL PD traits from both directories (14 files total)
# From munge_output: 8 files
# From munge_output_additional: 6 PD files
PD_FILES=(
    "${PD_DIR}/PD_females.A1ref.munged.sumstats.gz"
    "${PD_DIR}/PD_females.A1alt.munged.sumstats.gz"
    "${PD_DIR}/PD_males.A1ref.munged.sumstats.gz"
    "${PD_DIR}/PD_males.A1alt.munged.sumstats.gz"
    "${PD_DIR}/PDnp_females.A1ref.munged.sumstats.gz"
    "${PD_DIR}/PDnp_females.A1alt.munged.sumstats.gz"
    "${PD_DIR}/PDnp_males.A1ref.munged.sumstats.gz"
    "${PD_DIR}/PDnp_males.A1alt.munged.sumstats.gz"
    "${PD_ADDITIONAL_DIR}/PD_sexcombined_meta.tsv.A1ref.munged.sumstats.gz"
    "${PD_ADDITIONAL_DIR}/PD_sexcombined_meta.tsv.A1alt.munged.sumstats.gz"
    "${PD_ADDITIONAL_DIR}/NP_PD_sexcombined_meta.tsv.A1ref.munged.sumstats.gz"
    "${PD_ADDITIONAL_DIR}/NP_PD_sexcombined_meta.tsv.A1alt.munged.sumstats.gz"
    "${PD_ADDITIONAL_DIR}/PD_sexcomb.tsv.A1ref.munged.sumstats.gz"
    "${PD_ADDITIONAL_DIR}/PD_sexcomb.tsv.A1alt.munged.sumstats.gz"
    "${PD_ADDITIONAL_DIR}/PD_nalls2019_with_rsid.txt.A1ref.munged.sumstats.gz"
    "${PD_ADDITIONAL_DIR}/PD_nalls2019_with_rsid.txt.A1alt.munged.sumstats.gz"
)

PD_NAMES=(
    "PD_females.A1ref"
    "PD_females.A1alt"
    "PD_males.A1ref"
    "PD_males.A1alt"
    "PDnp_females.A1ref"
    "PDnp_females.A1alt"
    "PDnp_males.A1ref"
    "PDnp_males.A1alt"
    "PD_sexcombined_meta.A1ref"
    "PD_sexcombined_meta.A1alt"
    "NP_PD_sexcombined_meta.A1ref"
    "NP_PD_sexcombined_meta.A1alt"
    "PD_sexcomb.A1ref"
    "PD_sexcomb.A1alt"
    "PD_nalls2019.A1ref"
    "PD_nalls2019.A1alt"
)

# Define UACR traits (6 files: 3 phenotypes × 2 allele configs)
UACR_FILES=(
    "${UACR_DIR}/uacr_male.tsv.A1ref.munged.sumstats.gz"
    "${UACR_DIR}/uacr_male.tsv.A1alt.munged.sumstats.gz"
    "${UACR_DIR}/uacr_sex_combined.tsv.A1ref.munged.sumstats.gz"
    "${UACR_DIR}/uacr_sex_combined.tsv.A1alt.munged.sumstats.gz"
    "${UACR_DIR}/uacr_female.tsv.A1ref.munged.sumstats.gz"
    "${UACR_DIR}/uacr_female.tsv.A1alt.munged.sumstats.gz"
)

UACR_NAMES=(
    "uacr_male.A1ref"
    "uacr_male.A1alt"
    "uacr_sex_combined.A1ref"
    "uacr_sex_combined.A1alt"
    "uacr_female.A1ref"
    "uacr_female.A1alt"
)

# Calculate which PD and UACR combination this task should process
# Total combinations: 14 PD × 6 UACR = 84 combinations
# Each task processes one unique PD-UACR pair with their respective allele configurations

NUM_UACR=${#UACR_FILES[@]}
PD_INDEX=$((SLURM_ARRAY_TASK_ID / NUM_UACR))
UACR_INDEX=$((SLURM_ARRAY_TASK_ID % NUM_UACR))

PD_FILE="${PD_FILES[$PD_INDEX]}"
PD_NAME="${PD_NAMES[$PD_INDEX]}"
UACR_FILE="${UACR_FILES[$UACR_INDEX]}"
UACR_NAME="${UACR_NAMES[$UACR_INDEX]}"

echo "=========================================="
echo "Processing genetic correlation:"
echo "PD trait: ${PD_NAME}"
echo "UACR trait: ${UACR_NAME}"
echo "=========================================="
echo ""

# Check if input files exist
if [ ! -f "${PD_FILE}" ]; then
    echo "ERROR: PD file not found: ${PD_FILE}"
    exit 1
fi

if [ ! -f "${UACR_FILE}" ]; then
    echo "ERROR: UACR file not found: ${UACR_FILE}"
    exit 1
fi

# ============================================
# Run LDSC genetic correlation
# ============================================

OUTPUT_PREFIX="${OUTPUT_DIR}/rg_${PD_NAME}_vs_${UACR_NAME}"

echo "Running LDSC genetic correlation..."
echo "Started at: $(date)"
echo ""

time apptainer exec \
    --bind /home/lchang24:/home/lchang24 \
    --bind /scratch:/scratch \
    ${SIF_FILE} \
    ldsc.py \
    --rg ${PD_FILE},${UACR_FILE} \
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
echo "PD: ${PD_NAME}"
echo "UACR: ${UACR_NAME}"
echo "=========================================="