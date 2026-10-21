#!/bin/bash
#SBATCH --job-name=tax_769228
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G
#SBATCH --output=logs/tax_fixed_%j.out
#SBATCH --error=logs/tax_fixed_%j.err

module purge
module load bluebear
module load bear-apps/2022b
module load QIIME2/2024.2

cd /rds/projects/e/elhamsak-dissertation-autism/Microbiome/PRJNA769228 || exit 1
export TMPDIR=/rds/projects/e/elhamsak-dissertation-autism/tmp

CLASSIFIER=/rds/projects/e/elhamsak-dissertation-autism/Microbiome/silva-138-99-515-806-pretrained.qza

echo "=== Classify (mem-safe: n-jobs 1, batched): $(date) ==="
qiime feature-classifier classify-sklearn \
  --i-classifier $CLASSIFIER \
  --i-reads qiime2_results_fixed/rep-seqs.qza \
  --p-n-jobs 1 \
  --p-reads-per-batch 1000 \
  --o-classification qiime2_results_fixed/taxonomy.qza \
  --verbose

qiime metadata tabulate \
  --m-input-file qiime2_results_fixed/taxonomy.qza \
  --o-visualization qiime2_results_fixed/taxonomy.qzv

qiime taxa barplot \
  --i-table qiime2_results_fixed/table.qza \
  --i-taxonomy qiime2_results_fixed/taxonomy.qza \
  --m-metadata-file metadata/metadata_qiime2.tsv \
  --o-visualization qiime2_results_fixed/taxa-bar-plots.qzv

echo "=== Complete: $(date) ==="
