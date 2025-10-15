#!/bin/bash

#SBATCH --job-name=munge_array
#SBATCH --account=def-gsarah
#SBATCH --time=20:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --array=0-7
#SBATCH --output=munge_array-%A-%a.out
#SBATCH --error=munge_array-%A-%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=le.chang@umontreal.ca

# Exit on error and print commands
set -e
set -x

echo "=========================================="
echo "Array Job started at: $(date)"
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
SNP_LIST_PATH="/home/lchang24/projects/def-gsarah/lchang24/github/nf-genetic-correlations_fork/nf-genetic-correlations/data/ld_reference/w_hm3.snplist"
OUTPUT_DIR="/home/lchang24/scratch/munge_output"

# Create output directory if it doesn't exist
mkdir -p ${OUTPUT_DIR}

# Define all combinations: 4 files × 2 configurations = 8 array tasks
SUMSTATS_FILES=(
    "/home/lchang24/projects/def-gsarah/sumstats/kidney_neurodegen/PD_females.lava.gz"
    "/home/lchang24/projects/def-gsarah/sumstats/kidney_neurodegen/PD_females.lava.gz"
    "/home/lchang24/projects/def-gsarah/sumstats/kidney_neurodegen/PD_males.lava.gz"
    "/home/lchang24/projects/def-gsarah/sumstats/kidney_neurodegen/PD_males.lava.gz"
    "/home/lchang24/projects/def-gsarah/sumstats/kidney_neurodegen/PDnp_females.lava.gz"
    "/home/lchang24/projects/def-gsarah/sumstats/kidney_neurodegen/PDnp_females.lava.gz"
    "/home/lchang24/projects/def-gsarah/sumstats/kidney_neurodegen/PDnp_males.lava.gz"
    "/home/lchang24/projects/def-gsarah/sumstats/kidney_neurodegen/PDnp_males.lava.gz"
)

SAMPLE_SIZES=(
    104082
    104082
    110616
    110616
    19773
    19773
    24053
    24053
)

# Configuration: 0=A1ref, 1=A1alt (alternating pattern)
CONFIGS=(
    "A1ref"
    "A1alt"
    "A1ref"
    "A1alt"
    "A1ref"
    "A1alt"
    "A1ref"
    "A1alt"
)

# Get parameters for this array task
SUMSTATS="${SUMSTATS_FILES[$SLURM_ARRAY_TASK_ID]}"
N="${SAMPLE_SIZES[$SLURM_ARRAY_TASK_ID]}"
CONFIG="${CONFIGS[$SLURM_ARRAY_TASK_ID]}"

# Extract base filename
BASENAME=$(basename "$SUMSTATS" .lava.gz)

echo "=========================================="
echo "Processing: ${BASENAME}"
echo "Configuration: ${CONFIG}"
echo "Sample size: ${N}"
echo "=========================================="
echo ""

# ============================================
# Run munge_sumstats based on configuration
# ============================================

OUTPUT_PREFIX="${OUTPUT_DIR}/${BASENAME}.${CONFIG}.munged"

if [ "$CONFIG" == "A1ref" ]; then
    echo "Running with A1 as reference allele (--a1 A1 --a2 A2)"
    echo "Starting at: $(date)"
    echo ""
    
    time apptainer exec \
        --bind /home/lchang24:/home/lchang24 \
        --bind /scratch:/scratch \
        ${SIF_FILE} \
        munge_sumstats.py \
        --sumstats ${SUMSTATS} \
        --out ${OUTPUT_PREFIX} \
        --merge-alleles ${SNP_LIST_PATH} \
        --N ${N} \
        --ignore N \
        --snp SNP \
        --a1 A1 \
        --a2 A2 \
        --p P \
        --signed-sumstats BETA,0 \
        --chunksize 500000

else
    echo "Running with A1 as alternative allele (--a1 A2 --a2 A1)"
    echo "Starting at: $(date)"
    echo ""
    
    time apptainer exec \
        --bind /home/lchang24:/home/lchang24 \
        --bind /scratch:/scratch \
        ${SIF_FILE} \
        munge_sumstats.py \
        --sumstats ${SUMSTATS} \
        --out ${OUTPUT_PREFIX} \
        --merge-alleles ${SNP_LIST_PATH} \
        --N ${N} \
        --ignore N \
        --snp SNP \
        --a1 A2 \
        --a2 A1 \
        --p P \
        --signed-sumstats BETA,0 \
        --chunksize 500000
fi

# ============================================
# COMPLETION
# ============================================

echo ""
echo "================================================"
echo "COMPLETED!"
echo "================================================"
echo "Output: ${OUTPUT_PREFIX}.sumstats.gz"
echo "Log: ${OUTPUT_PREFIX}.log"
echo "Finished at: $(date)"
echo ""

# Show results
if [ -f "${OUTPUT_PREFIX}.sumstats.gz" ]; then
    echo "Output file:"
    ls -lh ${OUTPUT_PREFIX}.sumstats.gz
    echo ""
    echo "Number of SNPs:"
    zcat ${OUTPUT_PREFIX}.sumstats.gz | wc -l
else
    echo "Warning: Output file not found!"
    exit 1
fi

echo ""
echo "=========================================="
echo "Array task ${SLURM_ARRAY_TASK_ID} completed at: $(date)"
echo "=========================================="
