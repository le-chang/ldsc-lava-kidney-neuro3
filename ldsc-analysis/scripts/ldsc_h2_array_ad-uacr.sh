#!/bin/bash

#SBATCH --job-name=ldsc_h2_add
#SBATCH --account=def-gsarah
#SBATCH --time=20:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --array=0-13
#SBATCH --output=ldsc_h2_add-%A-%a.out
#SBATCH --error=ldsc_h2_add-%A-%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=le.chang@umontreal.ca

# Exit on error and print commands
set -e
set -x

echo "=========================================="
echo "LDSC h2 Job started at: $(date)"
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
MUNGED_DIR="/home/lchang24/scratch/munge_output_additional"
OUTPUT_DIR="/home/lchang24/scratch/ldsc_h2_results_additional"

# Create output directory if it doesn't exist
mkdir -p ${OUTPUT_DIR}

# Define munged sumstats files (14 files total)
MUNGED_FILES=(
    "${MUNGED_DIR}/PD_sexcombined_meta.tsv.A1ref.munged.sumstats.gz"
    "${MUNGED_DIR}/PD_sexcombined_meta.tsv.A1alt.munged.sumstats.gz"
    "${MUNGED_DIR}/NP_PD_sexcombined_meta.tsv.A1ref.munged.sumstats.gz"
    "${MUNGED_DIR}/NP_PD_sexcombined_meta.tsv.A1alt.munged.sumstats.gz"
    "${MUNGED_DIR}/PD_sexcomb.tsv.A1ref.munged.sumstats.gz"
    "${MUNGED_DIR}/PD_sexcomb.tsv.A1alt.munged.sumstats.gz"
    "${MUNGED_DIR}/PD_nalls2019_with_rsid.txt.A1ref.munged.sumstats.gz"
    "${MUNGED_DIR}/PD_nalls2019_with_rsid.txt.A1alt.munged.sumstats.gz"
    "${MUNGED_DIR}/uacr_male.tsv.A1ref.munged.sumstats.gz"
    "${MUNGED_DIR}/uacr_male.tsv.A1alt.munged.sumstats.gz"
    "${MUNGED_DIR}/uacr_sex_combined.tsv.A1ref.munged.sumstats.gz"
    "${MUNGED_DIR}/uacr_sex_combined.tsv.A1alt.munged.sumstats.gz"
    "${MUNGED_DIR}/uacr_female.tsv.A1ref.munged.sumstats.gz"
    "${MUNGED_DIR}/uacr_female.tsv.A1alt.munged.sumstats.gz"
)

# Get parameters for this array task
MUNGED_FILE="${MUNGED_FILES[$SLURM_ARRAY_TASK_ID]}"

# Extract base filename for output
BASENAME=$(basename "$MUNGED_FILE" .munged.sumstats.gz)

echo "=========================================="
echo "Processing: ${BASENAME}"
echo "Input file: ${MUNGED_FILE}"
echo "=========================================="
echo ""

# Check if input file exists
if [ ! -f "${MUNGED_FILE}" ]; then
    echo "ERROR: Input file not found: ${MUNGED_FILE}"
    exit 1
fi

# ============================================
# Run LDSC h2 estimation
# ============================================

OUTPUT_PREFIX="${OUTPUT_DIR}/${BASENAME}.h2"

echo "Running LDSC h2 estimation..."
echo "Started at: $(date)"
echo ""

time apptainer exec \
    --bind /home/lchang24:/home/lchang24 \
    --bind /scratch:/scratch \
    ${SIF_FILE} \
    ldsc.py \
    --h2 ${MUNGED_FILE} \
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
    echo "===== LDSC h2 Results ====="
    echo ""
    # Extract key results from log file
    grep -A 20 "Total Observed scale h2" ${OUTPUT_PREFIX}.log || echo "Results processing..."
    echo ""
    echo "Full log file: ${OUTPUT_PREFIX}.log"
else
    echo "Warning: Output log file not found!"
    exit 1
fi

echo ""
echo "=========================================="
echo "Array task ${SLURM_ARRAY_TASK_ID} completed at: $(date)"
echo "=========================================="
