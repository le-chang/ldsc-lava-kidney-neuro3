#!/bin/bash
#SBATCH --job-name=munge
#SBATCH --account=def-gsarah
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --array=0-12
#SBATCH --output=munge_%A_%a.out
#SBATCH --error=munge_%A_%a.err

# Simple parallel munge - based on what worked locally
# Just replicate the exact command that worked, for each trait

module load StdEnv/2023
module load apptainer/1.3.5

echo "=========================================="
echo "Array Job ${SLURM_ARRAY_TASK_ID}"
echo "Start: $(date)"
echo "=========================================="

SIF_FILE="/home/lchang24/projects/def-gsarah/lchang24/github/nf-genetic-correlations_fork/nf-genetic-correlations/bin/ldsc_latest.sif"
SNPLIST="/home/lchang24/projects/def-gsarah/lchang24/github/nf-genetic-correlations_fork/nf-genetic-correlations/data/ld_reference/w_hm3.snplist"
OUTPUT_DIR="/home/lchang24/scratch/ldsc_munged_all"

mkdir -p ${OUTPUT_DIR}

# Define traits: INPUT|OUTPUT|N|SNP|A1|A2|P|BETA
TRAITS=(
    "/home/lchang24/projects/def-gsarah/lchang24/github/nf-genetic-correlations_fork/nf-genetic-correlations/data/sumstats_original/AD_wightman2021_with_rsid.txt|AD_wightman2021|74004|SNP|A1|A2|P|BETA"
    "/home/lchang24/projects/def-gsarah/lchang24/github/nf-genetic-correlations_fork/nf-genetic-correlations/data/sumstats_original/AD_kunkle2019.txt|AD_kunkle2019|94437|MarkerName|Effect_allele|Non_Effect_allele|Pvalue|Beta"
    "/home/lchang24/projects/def-gsarah/lchang24/github/nf-genetic-correlations_fork/nf-genetic-correlations/data/sumstats_original/PD_nalls2019_with_rsid.txt|PD_nalls2019|482730|SNP|A1|A2|P|BETA"
    "/home/lchang24/projects/def-gsarah/sumstats/kidney_neurodegen/processed_sumstats/PD_sexcombined_meta.tsv|PD_sexcombined_meta|214698|SNP|A1|A2|P|BETA"
    "/home/lchang24/projects/def-gsarah/sumstats/kidney_neurodegen/processed_sumstats/PD_sexcomb.tsv|PD_sexcomb|463162|SNP|A1|A2|P|BETA"
    "/home/lchang24/projects/def-gsarah/sumstats/kidney_neurodegen/PD_males.lava.gz|PD_males|110616|SNP|A1|A2|P|BETA"
    "/home/lchang24/projects/def-gsarah/sumstats/kidney_neurodegen/PD_females.lava.gz|PD_females|104082|SNP|A1|A2|P|BETA"
    "/home/lchang24/projects/def-gsarah/sumstats/kidney_neurodegen/processed_sumstats/NP_PD_sexcombined_meta.tsv|NP_PD_sexcombined_meta|201289|SNP|A1|A2|P|BETA"
    "/home/lchang24/projects/def-gsarah/sumstats/kidney_neurodegen/PDnp_males.lava.gz|PDnp_males|24053|SNP|A1|A2|P|BETA"
    "/home/lchang24/projects/def-gsarah/sumstats/kidney_neurodegen/PDnp_females.lava.gz|PDnp_females|19773|SNP|A1|A2|P|BETA"
    "/home/lchang24/projects/def-gsarah/lchang24/github/nf-genetic-correlations_fork/nf-genetic-correlations/data/sumstats_original/uacr_sex_combined.tsv|uacr_sex_combined|122745|SNP|A1|A2|P|BETA"
    "/home/lchang24/projects/def-gsarah/lchang24/github/nf-genetic-correlations_fork/nf-genetic-correlations/data/sumstats_original/uacr_male.tsv|uacr_male|63642|SNP|A1|A2|P|BETA"
    "/home/lchang24/projects/def-gsarah/lchang24/github/nf-genetic-correlations_fork/nf-genetic-correlations/data/sumstats_original/uacr_female.tsv|uacr_female|59103|SNP|A1|A2|P|BETA"
)

# Get this task's trait
IFS='|' read -r INPUT OUTPUT N SNP A1 A2 P BETA <<< "${TRAITS[$SLURM_ARRAY_TASK_ID]}"

# Convert N to numeric (remove any non-numeric characters and ensure it's treated as number)
N=$(echo "$N" | sed 's/[^0-9]//g')

echo "Trait: ${OUTPUT}"
echo "N: ${N}"
echo ""

# A1A2 orientation - EXACT command that worked for you
echo "=== A1A2 ==="
apptainer run -W ${SLURM_TMPDIR} \
    ${SIF_FILE} munge_sumstats.py \
    --sumstats ${INPUT} \
    --merge-alleles ${SNPLIST} \
    --out ${OUTPUT_DIR}/${OUTPUT}_A1A2 \
    --N ${N} \
    --snp ${SNP} \
    --a1 ${A1} \
    --a2 ${A2} \
    --p ${P} \
    --signed-sumstats ${BETA},0

echo ""
echo "=== A2A1 ==="
apptainer run -W ${SLURM_TMPDIR} \
    ${SIF_FILE} munge_sumstats.py \
    --sumstats ${INPUT} \
    --merge-alleles ${SNPLIST} \
    --out ${OUTPUT_DIR}/${OUTPUT}_A2A1 \
    --N ${N} \
    --snp ${SNP} \
    --a1 ${A2} \
    --a2 ${A1} \
    --p ${P} \
    --signed-sumstats ${BETA},0

echo ""
echo "Done: $(date)"
