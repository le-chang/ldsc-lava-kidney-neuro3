#!/bin/bash

#SBATCH --time=2:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --job-name=diagnose_merge
#SBATCH --output=logs/diagnose_merge_%j.out
#SBATCH --error=logs/diagnose_merge_%j.err
#SBATCH --account=def-gsarah

set -x  # Print every command

echo "Starting diagnosis at: $(date)"

module load StdEnv/2023
module load apptainer/1.3.5

SIF_FILE="/home/lchang24/projects/def-gsarah/lchang24/github/nf-genetic-correlations_fork/nf-genetic-correlations/bin/ldsc_latest.sif"
OUTPUT_DIR="/home/lchang24/scratch/test_merge_diagnosis"
SNP_LIST_PATH="/home/lchang24/projects/def-gsarah/lchang24/github/nf-genetic-correlations_fork/nf-genetic-correlations/data/ld_reference/w_hm3.snplist"

mkdir -p ${OUTPUT_DIR}

# Use the smallest file first
SUMSTATS="/home/lchang24/projects/def-gsarah/sumstats/kidney_neurodegen/PDnp_females.lava.gz"

echo "================================================"
echo "Test 1: File checks"
echo "================================================"
echo "SUMSTATS file:"
ls -lh ${SUMSTATS}
echo ""
echo "SNP list file:"
ls -lh ${SNP_LIST_PATH}
wc -l ${SNP_LIST_PATH}
echo ""

echo "================================================"
echo "Test 2: WITHOUT --merge-alleles (should work fast)"
echo "================================================"
time apptainer exec \
    --bind /home/lchang24:/home/lchang24 \
    --bind /scratch:/scratch \
    ${SIF_FILE} \
    munge_sumstats.py \
    --sumstats ${SUMSTATS} \
    --out ${OUTPUT_DIR}/test_no_merge \
    --N 19773 \
    --ignore N \
    --snp SNP \
    --a1 A1 \
    --a2 A2 \
    --p P \
    --signed-sumstats BETA,0

echo "Test 2 completed at: $(date)"
echo ""

echo "================================================"
echo "Test 3: WITH --merge-alleles (this is where it might hang)"
echo "================================================"
echo "Starting at: $(date)"

# Add verbose output and use stdbuf to force unbuffered output
stdbuf -oL -eL apptainer exec \
    --bind /home/lchang24:/home/lchang24 \
    --bind /scratch:/scratch \
    ${SIF_FILE} \
    bash -c "
    set -x
    echo 'Inside container at: \$(date)'
    echo 'Checking files are accessible:'
    ls -lh ${SUMSTATS}
    ls -lh ${SNP_LIST_PATH}
    echo 'Starting munge_sumstats with merge-alleles at: \$(date)'
    
    time python /ldsc/munge_sumstats.py \
        --sumstats ${SUMSTATS} \
        --out ${OUTPUT_DIR}/test_with_merge \
        --merge-alleles ${SNP_LIST_PATH} \
        --N 19773 \
        --ignore N \
        --snp SNP \
        --a1 A1 \
        --a2 A2 \
        --p P \
        --signed-sumstats BETA,0
    
    echo 'Completed at: \$(date)'
    "

echo "================================================"
echo "Test 3 completed at: $(date)"
echo ""

echo "================================================"
echo "Test 4: Now try with LARGER file (PD_females)"
echo "================================================"
SUMSTATS_LARGE="/home/lchang24/projects/def-gsarah/sumstats/kidney_neurodegen/PD_females.lava.gz"
echo "File: ${SUMSTATS_LARGE}"
ls -lh ${SUMSTATS_LARGE}
echo "Starting at: $(date)"

time apptainer exec \
    --bind /home/lchang24:/home/lchang24 \
    --bind /scratch:/scratch \
    ${SIF_FILE} \
    munge_sumstats.py \
    --sumstats ${SUMSTATS_LARGE} \
    --out ${OUTPUT_DIR}/test_pd_females \
    --merge-alleles ${SNP_LIST_PATH} \
    --N 104082 \
    --ignore N \
    --snp SNP \
    --a1 A1 \
    --a2 A2 \
    --p P \
    --signed-sumstats BETA,0

echo "================================================"
echo "All tests completed at: $(date)"
echo "================================================"

echo "Output files created:"
ls -lh ${OUTPUT_DIR}/
