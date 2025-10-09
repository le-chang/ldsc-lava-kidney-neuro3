#!/bin/bash

#SBATCH --time=20:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --job-name=munge_sumstats
#SBATCH --output=logs/munge_%A_%a.out
#SBATCH --error=logs/munge_%A_%a.err
#SBATCH --account=def-gsarah

module load StdEnv/2020
module load python/2.7
module load scipy-stack/2020a

# LDSC munge_sumstats.py commands for all traits
# Running BOTH orientations (normal and flipped A1/A2) for each trait
# Delimiter: All files are tab-separated (\t)

#SIF_FILE="/home/lchang24/projects/def-gsarah/lchang24/github/nf-genetic-correlations_fork/nf-genetic-correlations/bin/ldsc_latest.sif"  # Update this path
OUTPUT_DIR="/home/lchang24/scratch/ldsc_munged_longer_local"
SNP_LIST_PATH="/home/lchang24/projects/def-gsarah/lchang24/github/nf-genetic-correlations_fork/nf-genetic-correlations/data/ld_reference/w_hm3.snplist"
LDSC_PATH="/home/lchang24/projects/def-gsarah/lchang24/software/ldsc/"
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
# AD_wightman2021 - Already done, but including both orientations
# ============================================================================
# Normal orientation
python ${LDSC_PATH}munge_sumstats.py \
    --sumstats /home/lchang24/projects/def-gsarah/lchang24/github/nf-genetic-correlations_fork/nf-genetic-correlations/data/sumstats_original/AD_wightman2021_with_rsid.txt \
    --out ${OUTPUT_DIR}/AD_wightman2021_normal \
    --merge-alleles ${SNP_LIST_PATH} \
    --N 762917 \
    --snp SNP \
    --a1 A1 \
    --a2 A2 \
    --p P \
    --signed-sumstats BETA,0

# Flipped orientation
python ${LDSC_PATH}munge_sumstats.py \
    --sumstats /home/lchang24/projects/def-gsarah/lchang24/github/nf-genetic-correlations_fork/nf-genetic-correlations/data/sumstats_original/AD_wightman2021_with_rsid.txt \
    --out ${OUTPUT_DIR}/AD_wightman2021_flipped \
    --merge-alleles ${SNP_LIST_PATH} \
    --N 762917 \
    --snp SNP \
    --a1 A2 \
    --a2 A1 \
    --p P \
    --signed-sumstats BETA,0

# ============================================================================
# PD_sexcombined_meta - Tab-delimited
# ============================================================================
# Normal orientation
python ${LDSC_PATH}munge_sumstats.py \
    --sumstats /home/lchang24/projects/def-gsarah/sumstats/kidney_neurodegen/processed_sumstats/PD_sexcombined_meta.tsv \
    --out ${OUTPUT_DIR}/PD_sexcombined_meta_normal \
    --merge-alleles ${SNP_LIST_PATH} \
    --N 214698 \
    --snp SNP \
    --a1 A1 \
    --a2 A2 \
    --p P \
    --signed-sumstats BETA,0

# Flipped orientation
python ${LDSC_PATH}munge_sumstats.py \
    --sumstats /home/lchang24/projects/def-gsarah/sumstats/kidney_neurodegen/processed_sumstats/PD_sexcombined_meta.tsv \
    --out ${OUTPUT_DIR}/PD_sexcombined_meta_flipped \
    --merge-alleles ${SNP_LIST_PATH} \
    --N 214698 \
    --snp SNP \
    --a1 A2 \
    --a2 A1 \
    --p P \
    --signed-sumstats BETA,0

# ============================================================================
# NP_PD_sexcombined_meta - Tab-delimited
# ============================================================================
# Normal orientation
python ${LDSC_PATH}munge_sumstats.py \
    --sumstats /home/lchang24/projects/def-gsarah/sumstats/kidney_neurodegen/processed_sumstats/NP_PD_sexcombined_meta.tsv \
    --out ${OUTPUT_DIR}/NP_PD_sexcombined_meta_normal \
    --merge-alleles ${SNP_LIST_PATH} \
    --N 201289 \
    --snp SNP \
    --a1 A1 \
    --a2 A2 \
    --p P \
    --signed-sumstats BETA,0

# Flipped orientation
python ${LDSC_PATH}munge_sumstats.py \
    --sumstats /home/lchang24/projects/def-gsarah/sumstats/kidney_neurodegen/processed_sumstats/NP_PD_sexcombined_meta.tsv \
    --out ${OUTPUT_DIR}/NP_PD_sexcombined_meta_flipped \
    --merge-alleles ${SNP_LIST_PATH} \
    --N 201289 \
    --snp SNP \
    --a1 A2 \
    --a2 A1 \
    --p P \
    --signed-sumstats BETA,0

# ============================================================================
# PD_sexcomb - Tab-delimited, has AF column
# ============================================================================
# Normal orientation
python ${LDSC_PATH}munge_sumstats.py \
    --sumstats /home/lchang24/projects/def-gsarah/sumstats/kidney_neurodegen/processed_sumstats/PD_sexcomb.tsv \
    --out ${OUTPUT_DIR}/PD_sexcomb_normal \
    --merge-alleles ${SNP_LIST_PATH} \
    --N 467813 \
    --snp SNP \
    --a1 A1 \
    --a2 A2 \
    --p P \
    --signed-sumstats BETA,0

# Flipped orientation
python ${LDSC_PATH}munge_sumstats.py \
    --sumstats /home/lchang24/projects/def-gsarah/sumstats/kidney_neurodegen/processed_sumstats/PD_sexcomb.tsv \
    --out ${OUTPUT_DIR}/PD_sexcomb_flipped \
    --merge-alleles ${SNP_LIST_PATH} \
    --N 467813 \
    --snp SNP \
    --a1 A2 \
    --a2 A1 \
    --p P \
    --signed-sumstats BETA,0

# ============================================================================
# PD_females.lava.gz - Gzipped, tab-delimited
# ============================================================================
# Normal orientation
python ${LDSC_PATH}munge_sumstats.py \
    --sumstats /home/lchang24/projects/def-gsarah/sumstats/kidney_neurodegen/PD_females.lava.gz \
    --out ${OUTPUT_DIR}/PD_females_normal \
    --merge-alleles ${SNP_LIST_PATH} \
    --N 104082 \
    --snp SNP \
    --a1 A1 \
    --a2 A2 \
    --p P \
    --signed-sumstats BETA,0

# Flipped orientation
python ${LDSC_PATH}munge_sumstats.py \
    --sumstats /home/lchang24/projects/def-gsarah/sumstats/kidney_neurodegen/PD_females.lava.gz \
    --out ${OUTPUT_DIR}/PD_females_flipped \
    --merge-alleles ${SNP_LIST_PATH} \
    --N 104082 \
    --snp SNP \
    --a1 A2 \
    --a2 A1 \
    --p P \
    --signed-sumstats BETA,0

# ============================================================================
# PD_males.lava.gz - Gzipped, tab-delimited
# ============================================================================
# Normal orientation
python ${LDSC_PATH}munge_sumstats.py \
    --sumstats /home/lchang24/projects/def-gsarah/sumstats/kidney_neurodegen/PD_males.lava.gz \
    --out ${OUTPUT_DIR}/PD_males_normal \
    --merge-alleles ${SNP_LIST_PATH} \
    --N 110616 \
    --snp SNP \
    --a1 A1 \
    --a2 A2 \
    --p P \
    --signed-sumstats BETA,0

# Flipped orientation
python ${LDSC_PATH}munge_sumstats.py \
    --sumstats /home/lchang24/projects/def-gsarah/sumstats/kidney_neurodegen/PD_males.lava.gz \
    --out ${OUTPUT_DIR}/PD_males_flipped \
    --merge-alleles ${SNP_LIST_PATH} \
    --N 110616 \
    --snp SNP \
    --a1 A2 \
    --a2 A1 \
    --p P \
    --signed-sumstats BETA,0

# ============================================================================
# PDnp_females.lava.gz - Gzipped, tab-delimited
# ============================================================================
# Normal orientation
python ${LDSC_PATH}munge_sumstats.py \
    --sumstats /home/lchang24/projects/def-gsarah/sumstats/kidney_neurodegen/PDnp_females.lava.gz \
    --out ${OUTPUT_DIR}/PDnp_females_normal \
    --merge-alleles ${SNP_LIST_PATH} \
    --N 19773 \
    --snp SNP \
    --a1 A1 \
    --a2 A2 \
    --p P \
    --signed-sumstats BETA,0

# Flipped orientation
python ${LDSC_PATH}munge_sumstats.py \
    --sumstats /home/lchang24/projects/def-gsarah/sumstats/kidney_neurodegen/PDnp_females.lava.gz \
    --out ${OUTPUT_DIR}/PDnp_females_flipped \
    --merge-alleles ${SNP_LIST_PATH} \
    --N 19773 \
    --snp SNP \
    --a1 A2 \
    --a2 A1 \
    --p P \
    --signed-sumstats BETA,0

# ============================================================================
# PDnp_males.lava.gz - Gzipped, tab-delimited
# ============================================================================
# Normal orientation
python ${LDSC_PATH}munge_sumstats.py \
    --sumstats /home/lchang24/projects/def-gsarah/sumstats/kidney_neurodegen/PDnp_males.lava.gz \
    --out ${OUTPUT_DIR}/PDnp_males_normal \
    --merge-alleles ${SNP_LIST_PATH} \
    --N 24053 \
    --snp SNP \
    --a1 A1 \
    --a2 A2 \
    --p P \
    --signed-sumstats BETA,0

# Flipped orientation
python ${LDSC_PATH}munge_sumstats.py \
    --sumstats /home/lchang24/projects/def-gsarah/sumstats/kidney_neurodegen/PDnp_males.lava.gz \
    --out ${OUTPUT_DIR}/PDnp_males_flipped \
    --merge-alleles ${SNP_LIST_PATH} \
    --N 24053 \
    --snp SNP \
    --a1 A2 \
    --a2 A1 \
    --p P \
    --signed-sumstats BETA,0

# ============================================================================
# PD_nalls2019_with_rsid.txt - Tab-delimited
# ============================================================================
# Normal orientation
python ${LDSC_PATH}munge_sumstats.py \
    --sumstats /home/lchang24/projects/def-gsarah/lchang24/github/nf-genetic-correlations_fork/nf-genetic-correlations/data/sumstats_original/PD_nalls2019_with_rsid.txt \
    --out ${OUTPUT_DIR}/PD_nalls2019_normal \
    --merge-alleles ${SNP_LIST_PATH} \
    --N 482730 \
    --snp SNP \
    --a1 A1 \
    --a2 A2 \
    --p P \
    --signed-sumstats BETA,0

# Flipped orientation
python ${LDSC_PATH}munge_sumstats.py \
    --sumstats /home/lchang24/projects/def-gsarah/lchang24/github/nf-genetic-correlations_fork/nf-genetic-correlations/data/sumstats_original/PD_nalls2019_with_rsid.txt \
    --out ${OUTPUT_DIR}/PD_nalls2019_flipped \
    --merge-alleles ${SNP_LIST_PATH} \
    --N 482730 \
    --snp SNP \
    --a1 A2 \
    --a2 A1 \
    --p P \
    --signed-sumstats BETA,0

# ============================================================================
# uacr_male.tsv - Tab-delimited, A2/A1 order reversed in file!
# ============================================================================
# Normal orientation (Note: file has A2 before A1 in header)
python ${LDSC_PATH}munge_sumstats.py \
    --sumstats /home/lchang24/projects/def-gsarah/lchang24/github/nf-genetic-correlations_fork/nf-genetic-correlations/data/sumstats_original/uacr_male.tsv \
    --out ${OUTPUT_DIR}/uacr_male_normal \
    --merge-alleles ${SNP_LIST_PATH} \
    --N 63642 \
    --snp SNP \
    --a1 A1 \
    --a2 A2 \
    --p P \
    --signed-sumstats BETA,0

# Flipped orientation
python ${LDSC_PATH}munge_sumstats.py \
    --sumstats /home/lchang24/projects/def-gsarah/lchang24/github/nf-genetic-correlations_fork/nf-genetic-correlations/data/sumstats_original/uacr_male.tsv \
    --out ${OUTPUT_DIR}/uacr_male_flipped \
    --merge-alleles ${SNP_LIST_PATH} \
    --N 63642 \
    --snp SNP \
    --a1 A2 \
    --a2 A1 \
    --p P \
    --signed-sumstats BETA,0

# ============================================================================
# uacr_sex_combined.tsv - Tab-delimited, A2/A1 order reversed in file!
# ============================================================================
# Normal orientation
python ${LDSC_PATH}munge_sumstats.py \
    --sumstats /home/lchang24/projects/def-gsarah/lchang24/github/nf-genetic-correlations_fork/nf-genetic-correlations/data/sumstats_original/uacr_sex_combined.tsv \
    --out ${OUTPUT_DIR}/uacr_sex_combined_normal \
    --merge-alleles ${SNP_LIST_PATH} \
    --N 122745 \
    --snp SNP \
    --a1 A1 \
    --a2 A2 \
    --p P \
    --signed-sumstats BETA,0

# Flipped orientation
python ${LDSC_PATH}munge_sumstats.py \
    --sumstats /home/lchang24/projects/def-gsarah/lchang24/github/nf-genetic-correlations_fork/nf-genetic-correlations/data/sumstats_original/uacr_sex_combined.tsv \
    --out ${OUTPUT_DIR}/uacr_sex_combined_flipped \
    --merge-alleles ${SNP_LIST_PATH} \
    --N 122745 \
    --snp SNP \
    --a1 A2 \
    --a2 A1 \
    --p P \
    --signed-sumstats BETA,0

# ============================================================================
# uacr_female.tsv - Tab-delimited, A2/A1 order reversed in file!
# ============================================================================
# Normal orientation
python ${LDSC_PATH}munge_sumstats.py \
    --sumstats /home/lchang24/projects/def-gsarah/lchang24/github/nf-genetic-correlations_fork/nf-genetic-correlations/data/sumstats_original/uacr_female.tsv \
    --out ${OUTPUT_DIR}/uacr_female_normal \
    --merge-alleles ${SNP_LIST_PATH} \
    --N 59103 \
    --snp SNP \
    --a1 A1 \
    --a2 A2 \
    --p P \
    --signed-sumstats BETA,0

# Flipped orientation
python ${LDSC_PATH}munge_sumstats.py \
    --sumstats /home/lchang24/projects/def-gsarah/lchang24/github/nf-genetic-correlations_fork/nf-genetic-correlations/data/sumstats_original/uacr_female.tsv \
    --out ${OUTPUT_DIR}/uacr_female_flipped \
    --merge-alleles ${SNP_LIST_PATH} \
    --N 59103 \
    --snp SNP \
    --a1 A2 \
    --a2 A1 \
    --p P \
    --signed-sumstats BETA,0

echo "All munge_sumstats.py commands completed!"
echo "Each trait has been processed in both normal and flipped A1/A2 orientations"
echo "Output files are in: ${OUTPUT_DIR}"