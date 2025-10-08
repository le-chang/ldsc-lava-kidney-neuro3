#!/bin/bash

#SBATCH --time=2:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --job-name=ldsc_rg_array
#SBATCH --array=1-108  # 6 UACR traits × (2 AD + 16 PD) = 108 total comparisons
#SBATCH --output=logs/ldsc_rg_%A_%a.out
#SBATCH --error=logs/ldsc_rg_%A_%a.err
#SBATCH --account=def-gsarah

module load StdEnv/2023
module load apptainer/1.3.5

# Create directories
mkdir -p logs
mkdir -p /home/lchang24/scratch/ldsc_results

# Set paths
SIF_FILE="/home/lchang24/projects/def-gsarah/lchang24/github/nf-genetic-correlations_fork/nf-genetic-correlations/bin/ldsc_latest.sif"
MUNGED_DIR="/home/lchang24/scratch/ldsc_munged"
OUTPUT_DIR="/home/lchang24/scratch/ldsc_results"

# LD Score reference panel - UPDATE THIS PATH!
LDSC_REF="/home/lchang24/projects/def-gsarah/lchang24/github/nf-genetic-correlations_fork/nf-genetic-correlations/data/ld_reference/eur_w_ld_chr/"  # UPDATE THIS PATH!
SNP_LIST_PATH="/home/lchang24/projects/def-gsarah/lchang24/github/nf-genetic-correlations_fork/nf-genetic-correlations/data/ld_reference/w_hm3.snplist"
# ============================================================================
# Define all trait pairs
# ============================================================================

# Create arrays for all comparisons
TRAIT1=()
TRAIT2=()

# UACR traits
UACR_TRAITS=(
    "uacr_male_normal"
    "uacr_male_flipped"
    "uacr_female_normal"
    "uacr_female_flipped"
    "uacr_sex_combined_normal"
    "uacr_sex_combined_flipped"
)

# AD traits
AD_TRAITS=(
    "AD_wightman2021_normal"
    "AD_wightman2021_flipped"
)

# PD traits
PD_TRAITS=(
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
)

# Build comparison pairs: UACR vs AD
for uacr in "${UACR_TRAITS[@]}"; do
    for ad in "${AD_TRAITS[@]}"; do
        TRAIT1+=("$uacr")
        TRAIT2+=("$ad")
    done
done

# Build comparison pairs: UACR vs PD
for uacr in "${UACR_TRAITS[@]}"; do
    for pd in "${PD_TRAITS[@]}"; do
        TRAIT1+=("$uacr")
        TRAIT2+=("$pd")
    done
done

# ============================================================================
# Run the correlation for this array task
# ============================================================================

# Get the trait pair for this array task (array indices start at 1)
IDX=$((SLURM_ARRAY_TASK_ID - 1))
TRAIT_A=${TRAIT1[$IDX]}
TRAIT_B=${TRAIT2[$IDX]}

echo "=========================================="
echo "Array Job ID: $SLURM_ARRAY_JOB_ID"
echo "Task ID: $SLURM_ARRAY_TASK_ID"
echo "Trait 1: $TRAIT_A"
echo "Trait 2: $TRAIT_B"
echo "Start Time: $(date)"
echo "=========================================="
echo ""

OUTPUT_NAME="${TRAIT_A}_vs_${TRAIT_B}"

# Run LDSC
apptainer run -W ${SLURM_TMPDIR} ${SIF_FILE} ldsc.py \
    --rg ${MUNGED_DIR}/${TRAIT_A}.sumstats.gz,${MUNGED_DIR}/${TRAIT_B}.sumstats.gz \
    --ref-ld-chr ${LDSC_REF} \
    --w-ld-chr ${LDSC_REF} \
    --out ${OUTPUT_DIR}/${OUTPUT_NAME}

echo ""
echo "=========================================="
echo "Completed: ${OUTPUT_NAME}"
echo "End Time: $(date)"
echo "=========================================="
