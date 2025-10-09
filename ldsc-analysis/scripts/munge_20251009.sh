#!/bin/bash

#SBATCH --time=2:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-24
#SBATCH --job-name=munge_array
#SBATCH --output=logs/munge_%A_%a.out
#SBATCH --error=logs/munge_%A_%a.err
#SBATCH --account=def-gsarah

set -e
set -o pipefail

echo "================================================"
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Job started at: $(date)"
echo "Running on node: $(hostname)"
echo "================================================"

module load StdEnv/2023
module load apptainer/1.3.5

SIF_FILE="/home/lchang24/projects/def-gsarah/lchang24/github/nf-genetic-correlations_fork/nf-genetic-correlations/bin/ldsc_latest.sif"
OUTPUT_DIR="/home/lchang24/scratch/ldsc_munged_longer"
SNP_LIST_PATH="/home/lchang24/projects/def-gsarah/lchang24/github/nf-genetic-correlations_fork/nf-genetic-correlations/data/ld_reference/w_hm3.snplist"

mkdir -p ${OUTPUT_DIR}

# Define all jobs in arrays
declare -a SUMSTATS_FILES=(
    "/home/lchang24/projects/def-gsarah/lchang24/github/nf-genetic-correlations_fork/nf-genetic-correlations/data/sumstats_original/AD_wightman2021_with_rsid.txt"
    "/home/lchang24/projects/def-gsarah/lchang24/github/nf-genetic-correlations_fork/nf-genetic-correlations/data/sumstats_original/AD_wightman2021_with_rsid.txt"
    "/home/lchang24/projects/def-gsarah/sumstats/kidney_neurodegen/processed_sumstats/PD_sexcombined_meta.tsv"
    "/home/lchang24/projects/def-gsarah/sumstats/kidney_neurodegen/processed_sumstats/PD_sexcombined_meta.tsv"
    "/home/lchang24/projects/def-gsarah/sumstats/kidney_neurodegen/processed_sumstats/NP_PD_sexcombined_meta.tsv"
    "/home/lchang24/projects/def-gsarah/sumstats/kidney_neurodegen/processed_sumstats/NP_PD_sexcombined_meta.tsv"
    "/home/lchang24/projects/def-gsarah/sumstats/kidney_neurodegen/processed_sumstats/PD_sexcomb.tsv"
    "/home/lchang24/projects/def-gsarah/sumstats/kidney_neurodegen/processed_sumstats/PD_sexcomb.tsv"
    "/home/lchang24/projects/def-gsarah/sumstats/kidney_neurodegen/PD_females.lava.gz"
    "/home/lchang24/projects/def-gsarah/sumstats/kidney_neurodegen/PD_females.lava.gz"
    "/home/lchang24/projects/def-gsarah/sumstats/kidney_neurodegen/PD_males.lava.gz"
    "/home/lchang24/projects/def-gsarah/sumstats/kidney_neurodegen/PD_males.lava.gz"
    "/home/lchang24/projects/def-gsarah/sumstats/kidney_neurodegen/PDnp_females.lava.gz"
    "/home/lchang24/projects/def-gsarah/sumstats/kidney_neurodegen/PDnp_females.lava.gz"
    "/home/lchang24/projects/def-gsarah/sumstats/kidney_neurodegen/PDnp_males.lava.gz"
    "/home/lchang24/projects/def-gsarah/sumstats/kidney_neurodegen/PDnp_males.lava.gz"
    "/home/lchang24/projects/def-gsarah/lchang24/github/nf-genetic-correlations_fork/nf-genetic-correlations/data/sumstats_original/PD_nalls2019_with_rsid.txt"
    "/home/lchang24/projects/def-gsarah/lchang24/github/nf-genetic-correlations_fork/nf-genetic-correlations/data/sumstats_original/PD_nalls2019_with_rsid.txt"
    "/home/lchang24/projects/def-gsarah/lchang24/github/nf-genetic-correlations_fork/nf-genetic-correlations/data/sumstats_original/uacr_male.tsv"
    "/home/lchang24/projects/def-gsarah/lchang24/github/nf-genetic-correlations_fork/nf-genetic-correlations/data/sumstats_original/uacr_male.tsv"
    "/home/lchang24/projects/def-gsarah/lchang24/github/nf-genetic-correlations_fork/nf-genetic-correlations/data/sumstats_original/uacr_sex_combined.tsv"
    "/home/lchang24/projects/def-gsarah/lchang24/github/nf-genetic-correlations_fork/nf-genetic-correlations/data/sumstats_original/uacr_sex_combined.tsv"
    "/home/lchang24/projects/def-gsarah/lchang24/github/nf-genetic-correlations_fork/nf-genetic-correlations/data/sumstats_original/uacr_female.tsv"
    "/home/lchang24/projects/def-gsarah/lchang24/github/nf-genetic-correlations_fork/nf-genetic-correlations/data/sumstats_original/uacr_female.tsv"
)

declare -a OUTPUT_NAMES=(
    "AD_wightman2021_normal"
    "AD_wightman2021_flipped"
    "PD_sexcombined_meta_normal"
    "PD_sexcombined_meta_flipped"
    "NP_PD_sexcombined_meta_normal"
    "NP_PD_sexcombined_meta_flipped"
    "PD_sexcomb_normal"
    "PD_sexcomb_flipped"
    "PD_females_normal"
    "PD_females_flipped"
    "PD_males_normal"
    "PD_males_flipped"
    "PDnp_females_normal"
    "PDnp_females_flipped"
    "PDnp_males_normal"
    "PDnp_males_flipped"
    "PD_nalls2019_normal"
    "PD_nalls2019_flipped"
    "uacr_male_normal"
    "uacr_male_flipped"
    "uacr_sex_combined_normal"
    "uacr_sex_combined_flipped"
    "uacr_female_normal"
    "uacr_female_flipped"
)

declare -a SAMPLE_SIZES=(
    762917 762917
    214698 214698
    201289 201289
    467813 467813
    104082 104082
    110616 110616
    19773 19773
    24053 24053
    482730 482730
    63642 63642
    122745 122745
    59103 59103
)

declare -a A1_COLS=(
    "A1" "A2"
    "A1" "A2"
    "A1" "A2"
    "A1" "A2"
    "A1" "A2"
    "A1" "A2"
    "A1" "A2"
    "A1" "A2"
    "A1" "A2"
    "A1" "A2"
    "A1" "A2"
    "A1" "A2"
)

declare -a A2_COLS=(
    "A2" "A1"
    "A2" "A1"
    "A2" "A1"
    "A2" "A1"
    "A2" "A1"
    "A2" "A1"
    "A2" "A1"
    "A2" "A1"
    "A2" "A1"
    "A2" "A1"
    "A2" "A1"
    "A2" "A1"
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

# Run munge_sumstats
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
echo "Output: ${OUTPUT_DIR}/${OUTPUT_NAME}.sumstats.gz"
ls -lh ${OUTPUT_DIR}/${OUTPUT_NAME}.sumstats.gz
echo "================================================"
