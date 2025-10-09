#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-8
#SBATCH --job-name=munge_lava
#SBATCH --output=logs/munge_lava_%A_%a.out
#SBATCH --error=logs/munge_lava_%A_%a.err
#SBATCH --account=def-gsarah

set -e
set -o pipefail
set -x  # ADD THIS - helps debug issues

echo "================================================"
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Job started at: $(date)"
echo "Running on node: $(hostname)"
echo "================================================"

module load StdEnv/2023
module load apptainer/1.3.5

SIF_FILE="/home/lchang24/projects/def-gsarah/lchang24/github/nf-genetic-correlations_fork/nf-genetic-correlations/bin/ldsc_latest.sif"
OUTPUT_DIR="/home/lchang24/scratch/ldsc_munged_longer_pd-2"
SNP_LIST_PATH="/home/lchang24/projects/def-gsarah/lchang24/github/nf-genetic-correlations_fork/nf-genetic-correlations/data/ld_reference/w_hm3.snplist"

mkdir -p ${OUTPUT_DIR}

# Define all jobs in arrays
declare -a SUMSTATS_FILES=(
    "/home/lchang24/projects/def-gsarah/sumstats/kidney_neurodegen/PD_females.lava.gz"
    "/home/lchang24/projects/def-gsarah/sumstats/kidney_neurodegen/PD_females.lava.gz"
    "/home/lchang24/projects/def-gsarah/sumstats/kidney_neurodegen/PD_males.lava.gz"
    "/home/lchang24/projects/def-gsarah/sumstats/kidney_neurodegen/PD_males.lava.gz"
    "/home/lchang24/projects/def-gsarah/sumstats/kidney_neurodegen/PDnp_females.lava.gz"
    "/home/lchang24/projects/def-gsarah/sumstats/kidney_neurodegen/PDnp_females.lava.gz"
    "/home/lchang24/projects/def-gsarah/sumstats/kidney_neurodegen/PDnp_males.lava.gz"
    "/home/lchang24/projects/def-gsarah/sumstats/kidney_neurodegen/PDnp_males.lava.gz"
)

declare -a OUTPUT_NAMES=(
    "PD_females_normal"
    "PD_females_flipped"
    "PD_males_normal"
    "PD_males_flipped"
    "PDnp_females_normal"
    "PDnp_females_flipped"
    "PDnp_males_normal"
    "PDnp_males_flipped"
)

declare -a SAMPLE_SIZES=(
    104082
    104082
    110616
    110616
    19773
    19773
    24053
    24053
)

declare -a A1_COLS=(
    "A1"
    "A2"
    "A1"
    "A2"
    "A1"
    "A2"
    "A1"
    "A2"
)

declare -a A2_COLS=(
    "A2"
    "A1"
    "A2"
    "A1"
    "A2"
    "A1"
    "A2"
    "A1"
)

# Get the parameters for this array task
IDX=$((SLURM_ARRAY_TASK_ID - 1))
SUMSTATS="${SUMSTATS_FILES[$IDX]}"
OUTPUT_NAME="${OUTPUT_NAMES[$IDX]}"
N="${SAMPLE_SIZES[$IDX]}"
A1="${A1_COLS[$IDX]}"
A2="${A2_COLS[$IDX]}"

echo "Processing file: ${SUMSTATS}"
echo "Output name: ${OUTPUT_NAME}"
echo "Sample size: ${N}"
echo "A1: ${A1}, A2: ${A2}"

# Check if input file exists
if [ ! -f "${SUMSTATS}" ]; then
    echo "ERROR: Input file not found: ${SUMSTATS}"
    exit 1
fi

echo "File exists, size: $(ls -lh ${SUMSTATS} | awk '{print $5}')"

# Peek at the file to see column names
echo "First few lines of input file:"
zcat ${SUMSTATS} | head -n 2

# Run munge_sumstats - simplified like the working diagnostic script
echo "Starting munge_sumstats at: $(date)"

apptainer exec ${SIF_FILE} munge_sumstats.py \
    --sumstats ${SUMSTATS} \
    --out ${OUTPUT_DIR}/${OUTPUT_NAME} \
    --merge-alleles ${SNP_LIST_PATH} \
    --N ${N} \
    --snp SNP \
    --a1 ${A1} \
    --a2 ${A2} \
    --p P \
    --signed-sumstats BETA,0

echo "================================================"
echo "Completed successfully at: $(date)"
if [ -f "${OUTPUT_DIR}/${OUTPUT_NAME}.sumstats.gz" ]; then
    echo "Output file created:"
    ls -lh ${OUTPUT_DIR}/${OUTPUT_NAME}.sumstats.gz
else
    echo "WARNING: Output file not found at ${OUTPUT_DIR}/${OUTPUT_NAME}.sumstats.gz"
    exit 1
fi
echo "================================================"
