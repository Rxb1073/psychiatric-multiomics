#!/bin/bash
#SBATCH --job-name=train_silva
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=/rds/projects/e/elhamsak-dissertation-autism/Microbiome/logs/train_silva_%j.out
#SBATCH --error=/rds/projects/e/elhamsak-dissertation-autism/Microbiome/logs/train_silva_%j.err

module purge
module load bear-apps/2023a
module load QIIME2/2025.4

export TMPDIR=/rds/projects/e/elhamsak-dissertation-autism/Microbiome/tmp
export TEMP=$TMPDIR
export TMP=$TMPDIR

mkdir -p /rds/projects/e/elhamsak-dissertation-autism/Microbiome/silva_training
mkdir -p /rds/projects/e/elhamsak-dissertation-autism/Microbiome/tmp
mkdir -p /rds/projects/e/elhamsak-dissertation-autism/Microbiome/logs

cd /rds/projects/e/elhamsak-dissertation-autism/Microbiome/silva_training

echo "=== SILVA Training Started: $(date) ==="

# Download SILVA 138 reference data
echo "Downloading SILVA sequences..."
wget -c https://data.qiime2.org/2022.8/common/silva-138-99-seqs.qza

echo "Downloading SILVA taxonomy..."
wget -c https://data.qiime2.org/2022.8/common/silva-138-99-tax.qza

# Extract V4 region (515F/806R)
echo "Extracting V4 region..."
qiime feature-classifier extract-reads \
  --i-sequences silva-138-99-seqs.qza \
  --p-f-primer GTGYCAGCMGCCGCGGTAA \
  --p-r-primer GGACTACNVGGGTWTCTAAT \
  --p-trunc-len 0 \
  --p-n-jobs 16 \
  --o-reads silva-138-99-515-806-seqs.qza

# Train classifier
echo "Training naive Bayes classifier..."
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads silva-138-99-515-806-seqs.qza \
  --i-reference-taxonomy silva-138-99-tax.qza \
  --o-classifier /rds/projects/e/elhamsak-dissertation-autism/Microbiome/silva-138-99-515-806-nb-classifier.qza

echo "=== SILVA Training Complete: $(date) ==="
ls -lh /rds/projects/e/elhamsak-dissertation-autism/Microbiome/silva-138-99-515-806-nb-classifier.qza