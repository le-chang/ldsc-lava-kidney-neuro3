#!/bin/bash

#SBATCH --job-name=munge_sumstats
#SBATCH --account=def-gsarah
#SBATCH --time=20:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=le.chang@umontreal.ca

# Exit on error and print commands
set -e
set -x

echo "=========================================="
echo "Job started at: $(date)"
echo "Job ID: ${SLURM_JOB_ID}"
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
FILTERED_SUMSTATS="/home/lchang24/scratch/test_merge_diagnosis/PDnp_females.lava.filtered.gz"
OUTPUT_DIR="/home/lchang24/scratch/test_merge_diagnosis"
N=19773

# Create output directory if it doesn't exist
mkdir -p ${OUTPUT_DIR}

# ============================================
# Run munge_sumstats with --merge-alleles
# ============================================

echo "Starting munge with --merge-alleles at: $(date)"
echo "Input: ${FILTERED_SUMSTATS}"
echo "Reference: ${SNP_LIST_PATH}"
echo "Sample size: ${N}"
echo ""

OUTPUT_PREFIX="${OUTPUT_DIR}/PDnp_females.munged"

echo "Running munge_sumstats..."
time apptainer exec \
    --bind /home/lchang24:/home/lchang24 \
    --bind /scratch:/scratch \
    ${SIF_FILE} \
    munge_sumstats.py \
    --sumstats ${FILTERED_SUMSTATS} \
    --out ${OUTPUT_PREFIX} \
    --merge-alleles ${SNP_LIST_PATH} \
    --N ${N} \
    --ignore N \
    --snp SNP \
    --a1 A1 \
    --a2 A2 \
    --p P \
    --signed-sumstats BETA,0

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
echo "Job completed at: $(date)"
echo "=========================================="
