#!/bin/bash

#SBATCH --job-name=ldsc_rg_pdnp_sc
#SBATCH --account=def-gsarah
#SBATCH --time=3:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --array=0-11
#SBATCH --output=ldsc_rg_pdnp_sc-%A-%a.out
#SBATCH --error=ldsc_rg_pdnp_sc-%A-%a.err
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
PDNP_DIR="/home/lchang24/scratch/munge_output_additional"
UACR_DIR="/home/lchang24/scratch/munge_output_additional"
OUTPUT_DIR="/home/lchang24/scratch/ldsc_rg_results_pdnp_sexcomb"

# Create output directory if it doesn't exist
mkdir -p ${OUTPUT_DIR}

# Define PDnp_sexcombined traits (2 files: 1 phenotype × 2 allele configs)
PDNP_FILES=(
    "${PDNP_DIR}/PDnp_sexcombined.A1ref.munged.sumstats.gz"
    "${PDNP_DIR}/PDnp_sexcombined.A1alt.munged.sumstats.gz"
    "${PDNP_DIR}/PDnp_sexcombined.A1ref.munged.sumstats.gz"
    "${PDNP_DIR}/PDnp_sexcombined.A1alt.munged.sumstats.gz"
    "${PDNP_DIR}/PDnp_sexcombined.A1ref.munged.sumstats.gz"
    "${PDNP_DIR}/PDnp_sexcombined.A1alt.munged.sumstats.gz"
    "${PDNP_DIR}/PDnp_sexcombined.A1ref.munged.sumstats.gz"
    "${PDNP_DIR}/PDnp_sexcombined.A1alt.munged.sumstats.gz"
    "${PDNP_DIR}/PDnp_sexcombined.A1ref.munged.sumstats.gz"
    "${PDNP_DIR}/PDnp_sexcombined.A1alt.munged.sumstats.gz"
    "${PDNP_DIR}/PDnp_sexcombined.A1ref.munged.sumstats.gz"
    "${PDNP_DIR}/PDnp_sexcombined.A1alt.munged.sumstats.gz"
)

PDNP_NAMES=(
    "PDnp_sexcombined.A1ref"
    "PDnp_sexcombined.A1alt"
    "PDnp_sexcombined.A1ref"
    "PDnp_sexcombined.A1alt"
    "PDnp_sexcombined.A1ref"
    "PDnp_sexcombined.A1alt"
    "PDnp_sexcombined.A1ref"
    "PDnp_sexcombined.A1alt"
    "PDnp_sexcombined.A1ref"
    "PDnp_sexcombined.A1alt"
    "PDnp_sexcombined.A1ref"
    "PDnp_sexcombined.A1alt"
)

# Define UACR traits (6 files: 3 phenotypes × 2 allele configs)
UACR_FILES=(
    "${UACR_DIR}/uacr_male.tsv.A1ref.munged.sumstats.gz"
    "${UACR_DIR}/uacr_male.tsv.A1alt.munged.sumstats.gz"
    "${UACR_DIR}/uacr_sex_combined.tsv.A1ref.munged.sumstats.gz"
    "${UACR_DIR}/uacr_sex_combined.tsv.A1alt.munged.sumstats.gz"
    "${UACR_DIR}/uacr_female.tsv.A1ref.munged.sumstats.gz"
    "${UACR_DIR}/uacr_female.tsv.A1alt.munged.sumstats.gz"
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
    "uacr_male.A1ref"
    "uacr_male.A1alt"
    "uacr_sex_combined.A1ref"
    "uacr_sex_combined.A1alt"
    "uacr_female.A1ref"
    "uacr_female.A1alt"
)

# Get the specific combination for this task
PDNP_FILE="${PDNP_FILES[$SLURM_ARRAY_TASK_ID]}"
PDNP_NAME="${PDNP_NAMES[$SLURM_ARRAY_TASK_ID]}"
UACR_FILE="${UACR_FILES[$SLURM_ARRAY_TASK_ID]}"
UACR_NAME="${UACR_NAMES[$SLURM_ARRAY_TASK_ID]}"

echo "=========================================="
echo "Processing genetic correlation:"
echo "PDnp trait: ${PDNP_NAME}"
echo "UACR trait: ${UACR_NAME}"
echo "=========================================="
echo ""

# Check if input files exist
if [ ! -f "${PDNP_FILE}" ]; then
    echo "ERROR: PDnp file not found: ${PDNP_FILE}"
    exit 1
fi

if [ ! -f "${UACR_FILE}" ]; then
    echo "ERROR: UACR file not found: ${UACR_FILE}"
    exit 1
fi

# ============================================
# Run LDSC genetic correlation
# ============================================

OUTPUT_PREFIX="${OUTPUT_DIR}/rg_${PDNP_NAME}_vs_${UACR_NAME}"

echo "Running LDSC genetic correlation..."
echo "Started at: $(date)"
echo ""

time apptainer exec \
    --bind /home/lchang24:/home/lchang24 \
    --bind /scratch:/scratch \
    ${SIF_FILE} \
    ldsc.py \
    --rg ${PDNP_FILE},${UACR_FILE} \
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
echo "PDnp: ${PDNP_NAME}"
echo "UACR: ${UACR_NAME}"
echo "=========================================="
