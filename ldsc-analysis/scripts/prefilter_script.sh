#!/bin/bash

set -e  # Exit on error
set -x  # Print commands

echo "Starting pre-filtering at: $(date)"

# ============================================
# CONFIGURATION
# ============================================
SNP_LIST_PATH="/home/lchang24/projects/def-gsarah/lchang24/github/nf-genetic-correlations_fork/nf-genetic-correlations/data/ld_reference/w_hm3.snplist"
SUMSTATS="/home/lchang24/projects/def-gsarah/sumstats/kidney_neurodegen/PDnp_females.lava.gz"
OUTPUT_DIR="/home/lchang24/scratch/test_merge_diagnosis"

mkdir -p ${OUTPUT_DIR}

# ============================================
# Pre-filter: Keep only HapMap3 SNPs
# ============================================
echo "Input file: ${SUMSTATS}"
echo "Reference: ${SNP_LIST_PATH}"
echo ""

FILTERED_OUTPUT="${OUTPUT_DIR}/PDnp_females.lava.filtered.gz"

echo "Filtering to HapMap3 SNPs..."

zcat ${SUMSTATS} | awk '
BEGIN {
    # Load all HapMap3 SNP IDs into memory
    while ((getline < "'${SNP_LIST_PATH}'") > 0) {
        snps[$1] = 1
    }
    close("'${SNP_LIST_PATH}'")
    print "Loaded " length(snps) " reference SNPs" > "/dev/stderr"
}
NR == 1 {
    # Keep header
    print $0
    next
}
{
    # Keep SNP if it exists in reference
    if ($1 in snps) {
        print $0
    }
}' | gzip > ${FILTERED_OUTPUT}

echo ""
echo "Done!"
echo "Original SNPs: $(zcat ${SUMSTATS} | wc -l)"
echo "Filtered SNPs: $(zcat ${FILTERED_OUTPUT} | wc -l)"
echo "Output: ${FILTERED_OUTPUT}"
echo "Finished at: $(date)"
