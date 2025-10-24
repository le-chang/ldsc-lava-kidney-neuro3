#!/bin/bash

#SBATCH --job-name=munge_pdnp_sc
#SBATCH --account=def-gsarah
#SBATCH --time=8:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --array=0-1
#SBATCH --output=munge_pdnp_sc-%A-%a.out
#SBATCH --error=munge_pdnp_sc-%A-%a.err
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
OUTPUT_DIR="/home/lchang24/scratch/munge_output_additional"

# Create output directory if it doesn't exist
mkdir -p ${OUTPUT_DIR}

# Define the sumstats file and sample size
SUMSTATS="/home/lchang24/projects/def-gsarah/lchang24/github/nf-genetic-correlations_fork/nf-genetic-correlations/data/sumstats_original/PDnp_sexcombined.tsv.gz"
N=1437688
BASENAME="PDnp_sexcombined"

# Configuration: 0=A1ref, 1=A1alt
CONFIGS=(
    "A1ref"
    "A1alt"
)

# Get configuration for this array task
CONFIG="${CONFIGS[$SLURM_ARRAY_TASK_ID]}"

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
    echo "Running with effect_allele as reference allele (--a1 effect_allele --a2 other_allele)"
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
        --snp rsid \
        --a1 effect_allele \
        --a2 other_allele \
        --p p_value \
        --signed-sumstats beta,0 \
        --chunksize 500000

else
    echo "Running with effect_allele as alternative allele (--a1 other_allele --a2 effect_allele)"
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
        --snp rsid \
        --a1 other_allele \
        --a2 effect_allele \
        --p p_value \
        --signed-sumstats beta,0 \
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
