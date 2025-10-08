#!/bin/bash

#SBATCH --time=20:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --job-name=munge_lava_fix
#SBATCH --output=logs/munge_lava_%A_%a.out
#SBATCH --error=logs/munge_lava_%A_%a.err
#SBATCH --account=def-gsarah

module load StdEnv/2020
module load python/2.7
module load scipy-stack/2020a
#SIF_FILE="/home/lchang24/projects/def-gsarah/lchang24/github/nf-genetic-correlations_fork/nf-genetic-correlations/bin/ldsc_latest.sif"
OUTPUT_DIR="/home/lchang24/scratch/ldsc_munged_longer_local"
SNP_LIST_PATH="/home/lchang24/projects/def-gsarah/lchang24/github/nf-genetic-correlations_fork/nf-genetic-correlations/data/ld_reference/w_hm3.snplist"
LDSC_PATH="/home/lchang24/projects/def-gsarah/lchang24/software/ldsc"
mkdir -p ${OUTPUT_DIR}

virtualenv --no-download $SLURM_TMPDIR/env
source $SLURM_TMPDIR/env/bin/activate
pip install --upgrade pip
# pip install bitarray
# pip install pybedtools==0.7.10
pip install contextlib2==0.2
pip install importlib_metadata==1.5.0
pip install jsonschema==2.6.0
pip install six==1.16.0+computecanada
pip install bitarray==0.8.0

pip install -r ${LDSC_PATH}/ldsc_requirements.txt


# ============================================================================
# FIX: Add --ignore N to tell LDSC to ignore the N column from the file
# ============================================================================

# ============================================================================
# PD_females.lava.gz - FIXED
# ============================================================================
echo "Processing PD_females - normal orientation..."
python ${LDSC_PATH}munge_sumstats.py \
    --sumstats /home/lchang24/projects/def-gsarah/sumstats/kidney_neurodegen/PD_females.lava.gz \
    --out ${OUTPUT_DIR}/PD_females_normal \
    --merge-alleles ${SNP_LIST_PATH} \
    --N 104082 \
    --ignore N \
    --snp SNP \
    --a1 A1 \
    --a2 A2 \
    --p P \
    --signed-sumstats BETA,0

echo "Processing PD_females - flipped orientation..."
python ${LDSC_PATH}munge_sumstats.py \
    --sumstats /home/lchang24/projects/def-gsarah/sumstats/kidney_neurodegen/PD_females.lava.gz \
    --out ${OUTPUT_DIR}/PD_females_flipped \
    --merge-alleles ${SNP_LIST_PATH} \
    --N 104082 \
    --ignore N \
    --snp SNP \
    --a1 A2 \
    --a2 A1 \
    --p P \
    --signed-sumstats BETA,0

# ============================================================================
# PD_males.lava.gz - FIXED
# ============================================================================
echo "Processing PD_males - normal orientation..."
python ${LDSC_PATH}munge_sumstats.py \
    --sumstats /home/lchang24/projects/def-gsarah/sumstats/kidney_neurodegen/PD_males.lava.gz \
    --out ${OUTPUT_DIR}/PD_males_normal \
    --merge-alleles ${SNP_LIST_PATH} \
    --N 110616 \
    --ignore N \
    --snp SNP \
    --a1 A1 \
    --a2 A2 \
    --p P \
    --signed-sumstats BETA,0

echo "Processing PD_males - flipped orientation..."
python ${LDSC_PATH}munge_sumstats.py \
    --sumstats /home/lchang24/projects/def-gsarah/sumstats/kidney_neurodegen/PD_males.lava.gz \
    --out ${OUTPUT_DIR}/PD_males_flipped \
    --merge-alleles ${SNP_LIST_PATH} \
    --N 110616 \
    --ignore N \
    --snp SNP \
    --a1 A2 \
    --a2 A1 \
    --p P \
    --signed-sumstats BETA,0

# ============================================================================
# PDnp_females.lava.gz - FIXED
# ============================================================================
echo "Processing PDnp_females - normal orientation..."
python ${LDSC_PATH}munge_sumstats.py \
    --sumstats /home/lchang24/projects/def-gsarah/sumstats/kidney_neurodegen/PDnp_females.lava.gz \
    --out ${OUTPUT_DIR}/PDnp_females_normal \
    --merge-alleles ${SNP_LIST_PATH} \
    --N 19773 \
    --ignore N \
    --snp SNP \
    --a1 A1 \
    --a2 A2 \
    --p P \
    --signed-sumstats BETA,0

echo "Processing PDnp_females - flipped orientation..."
python ${LDSC_PATH}munge_sumstats.py \
    --sumstats /home/lchang24/projects/def-gsarah/sumstats/kidney_neurodegen/PDnp_females.lava.gz \
    --out ${OUTPUT_DIR}/PDnp_females_flipped \
    --merge-alleles ${SNP_LIST_PATH} \
    --N 19773 \
    --ignore N \
    --snp SNP \
    --a1 A2 \
    --a2 A1 \
    --p P \
    --signed-sumstats BETA,0

# ============================================================================
# PDnp_males.lava.gz - FIXED
# ============================================================================
echo "Processing PDnp_males - normal orientation..."
python ${LDSC_PATH}munge_sumstats.py \
    --sumstats /home/lchang24/projects/def-gsarah/sumstats/kidney_neurodegen/PDnp_males.lava.gz \
    --out ${OUTPUT_DIR}/PDnp_males_normal \
    --merge-alleles ${SNP_LIST_PATH} \
    --N 24053 \
    --ignore N \
    --snp SNP \
    --a1 A1 \
    --a2 A2 \
    --p P \
    --signed-sumstats BETA,0

echo "Processing PDnp_males - flipped orientation..."
python ${LDSC_PATH}munge_sumstats.py \
    --sumstats /home/lchang24/projects/def-gsarah/sumstats/kidney_neurodegen/PDnp_males.lava.gz \
    --out ${OUTPUT_DIR}/PDnp_males_flipped \
    --merge-alleles ${SNP_LIST_PATH} \
    --N 24053 \
    --ignore N \
    --snp SNP \
    --a1 A2 \
    --a2 A1 \
    --p P \
    --signed-sumstats BETA,0

echo ""
echo "=========================================="
echo "All .lava.gz files processed successfully!"
echo "End Time: $(date)"
echo "=========================================="
