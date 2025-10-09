#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-8
#SBATCH --job-name=munge_lava
#SBATCH --output=logs/munge_lava_%A_%a.out
#SBATCH --error=logs/munge_lava_%A_%a.err
#SBATCH --account=def-gsarah

export APPTAINER_BIND="/home/lchang24,/scratch,${SLURM_TMPDIR}"

echo "================================================"
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Job started at: $(date)"
echo "================================================"

module load StdEnv/2023
module load apptainer/1.3.5

SIF_FILE="/home/lchang24/projects/def-gsarah/lchang24/github/nf-genetic-correlations_fork/nf-genetic-correlations/bin/ldsc_latest.sif"
OUTPUT_DIR="/home/lchang24/scratch/ldsc_munged_longer2"
SNP_LIST_PATH="/home/lchang24/projects/def-gsarah/lchang24/github/nf-genetic-correlations_fork/nf-genetic-correlations/data/ld_reference/w_hm3.snplist"

mkdir -p ${OUTPUT_DIR}

# Arrays as before...
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
    "PD_females_normal" "PD_females_flipped"
    "PD_males_normal" "PD_males_flipped"
    "PDnp_females_normal" "PDnp_females_flipped"
    "PDnp_males_normal" "PDnp_males_flipped"
)

declare -a SAMPLE_SIZES=(104082 104082 110616 110616 19773 19773 24053 24053)
declare -a A1_COLS=("A1" "A2" "A1" "A2" "A1" "A2" "A1" "A2")
declare -a A2_COLS=("A2" "A1" "A2" "A1" "A2" "A1" "A2" "A1")

IDX=$((SLURM_ARRAY_TASK_ID - 1))
SUMSTATS="${SUMSTATS_FILES[$IDX]}"
OUTPUT_NAME="${OUTPUT_NAMES[$IDX]}"
N="${SAMPLE_SIZES[$IDX]}"
A1="${A1_COLS[$IDX]}"
A2="${A2_COLS[$IDX]}"

echo "Processing: ${OUTPUT_NAME}"
echo "Started at: $(date)"

apptainer exec ${SIF_FILE} munge_sumstats.py \
    --sumstats ${SUMSTATS} \
    --out ${OUTPUT_DIR}/${OUTPUT_NAME} \
    --merge-alleles ${SNP_LIST_PATH} \
    --N ${N} \
    --ignore N \
    --snp SNP \
    --a1 ${A1} \
    --a2 ${A2} \
    --p P \
    --signed-sumstats BETA,0

echo "Completed at: $(date)"
ls -lh ${OUTPUT_DIR}/${OUTPUT_NAME}.sumstats.gz
