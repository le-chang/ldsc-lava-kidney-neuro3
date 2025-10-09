#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --job-name=munge_lava_single
#SBATCH --output=logs/munge_lava_single.out
#SBATCH --error=logs/munge_lava_single.err
#SBATCH --account=def-gsarah

set -e
set -o pipefail
set -x  # Print each command for debugging

echo "================================================"
echo "Job started at: $(date)"
echo "Running on node: $(hostname)"
echo "================================================"

module load StdEnv/2023
module load apptainer/1.3.5

SIF_FILE="/home/lchang24/projects/def-gsarah/lchang24/github/nf-genetic-correlations_fork/nf-genetic-correlations/bin/ldsc_latest.sif"
OUTPUT_DIR="/home/lchang24/scratch/ldsc_munged_longer_pd334"
SNP_LIST_PATH="/home/lchang24/projects/def-gsarah/lchang24/github/nf-genetic-correlations_fork/nf-genetic-correlations/data/ld_reference/w_hm3.snplist"

mkdir -p ${OUTPUT_DIR}
mkdir -p logs

# Single job example - PD_females_normal
SUMSTATS="/home/lchang24/projects/def-gsarah/sumstats/kidney_neurodegen/PD_females.lava.gz"
OUTPUT_NAME="PD_females_normal"
N=104082
A1="A1"
A2="A2"

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

# Check if SNP list exists
if [ ! -f "${SNP_LIST_PATH}" ]; then
    echo "ERROR: SNP list not found: ${SNP_LIST_PATH}"
    exit 1
fi

echo "SNP list exists, size: $(ls -lh ${SNP_LIST_PATH} | awk '{print $5}')"

# Run munge_sumstats - following the working diagnostic script pattern
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
