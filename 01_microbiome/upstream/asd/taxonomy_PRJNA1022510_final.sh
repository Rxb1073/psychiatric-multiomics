#!/bin/bash
#SBATCH --job-name=tax_PRJNA1022510
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=logs/tax_final_%j.out
#SBATCH --error=logs/tax_final_%j.err

module purge
module load bluebear
module load bear-apps/2022b
module load QIIME2/2024.2

cd /rds/projects/e/elhamsak-dissertation-autism/Microbiome/PRJNA1022510 || exit 1
export TMPDIR=/rds/projects/e/elhamsak-dissertation-autism/tmp

CLASSIFIER=/rds/projects/e/elhamsak-dissertation-autism/Microbiome/silva-138-99-515-806-pretrained.qza

# --- 1. Promote final files to canonical names ---
cp qiime2_results_sra/table_final.qza      qiime2_results_sra/table.qza
cp qiime2_results_sra/rep-seqs_final.qza   qiime2_results_sra/rep-seqs.qza

echo "=== Taxonomic classification (SILVA 138, consistent with other 3 datasets): $(date) ==="
qiime feature-classifier classify-sklearn \
  --i-classifier $CLASSIFIER \
  --i-reads qiime2_results_sra/rep-seqs.qza \
  --p-n-jobs 16 \
  --o-classification qiime2_results_sra/taxonomy.qza \
  --verbose

qiime metadata tabulate \
  --m-input-file qiime2_results_sra/taxonomy.qza \
  --o-visualization qiime2_results_sra/taxonomy.qzv

echo "=== Taxa bar plots: $(date) ==="
qiime taxa barplot \
  --i-table qiime2_results_sra/table.qza \
  --i-taxonomy qiime2_results_sra/taxonomy.qza \
  --m-metadata-file metadata/GSE244391_metadata.tsv \
  --o-visualization qiime2_results_sra/taxa-bar-plots.qzv

echo "=== Feature table summary: $(date) ==="
qiime feature-table summarize \
  --i-table qiime2_results_sra/table.qza \
  --m-sample-metadata-file metadata/GSE244391_metadata.tsv \
  --o-visualization qiime2_results_sra/table-summary.qzv

echo "=== Complete: $(date) ==="