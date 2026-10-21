#!/bin/bash
#SBATCH --job-name=taxonomy_PRJNA1022510
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --output=logs/taxonomy_%j.out
#SBATCH --error=logs/taxonomy_%j.err

module purge
module load bear-apps/2022b/live
module load QIIME2/2024.2

export TMPDIR=/rds/projects/e/elhamsak-dissertation-autism/tmp
mkdir -p $TMPDIR

cd /rds/projects/e/elhamsak-dissertation-autism/Microbiome/PRJNA1022510 || exit 1

CLASSIFIER=/rds/projects/e/elhamsak-dissertation-autism/Microbiome/silva-138-99-515-806-pretrained.qza

echo "=== Taxonomy Started: $(date) ==="

qiime feature-classifier classify-sklearn \
  --i-classifier $CLASSIFIER \
  --i-reads qiime2_results/rep-seqs.qza \
  --p-n-jobs 16 \
  --o-classification qiime2_results/taxonomy.qza

qiime metadata tabulate \
  --m-input-file qiime2_results/taxonomy.qza \
  --o-visualization qiime2_results/taxonomy.qzv

qiime taxa barplot \
  --i-table qiime2_results/table.qza \
  --i-taxonomy qiime2_results/taxonomy.qza \
  --m-metadata-file /rds/projects/e/elhamsak-dissertation-autism/Microbiome/PRJNA1022510/metadata/metadata_qiime2.tsv \
  --o-visualization qiime2_results/taxa-bar-plots.qzv

echo "=== Taxonomy Complete: $(date) ==="