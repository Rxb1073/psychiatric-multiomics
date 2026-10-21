#!/bin/bash
#SBATCH --job-name=tax_1022_fix
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G
#SBATCH --output=logs/tax_oomfix_%j.out
#SBATCH --error=logs/tax_oomfix_%j.err

module purge
module load bluebear
module load bear-apps/2022b
module load QIIME2/2024.2

cd /rds/projects/e/elhamsak-dissertation-autism/Microbiome/PRJNA1022510 || exit 1
export TMPDIR=/rds/projects/e/elhamsak-dissertation-autism/tmp

CLASSIFIER=/rds/projects/e/elhamsak-dissertation-autism/Microbiome/silva-138-99-515-806-pretrained.qza

echo "=== Classify (mem=120G, n-jobs=1 to limit memory): $(date) ==="
qiime feature-classifier classify-sklearn \
  --i-classifier $CLASSIFIER \
  --i-reads qiime2_results_sra/rep-seqs_final.qza \
  --p-n-jobs 1 \
  --p-reads-per-batch 1000 \
  --o-classification qiime2_results_sra/taxonomy.qza \
  --verbose

qiime metadata tabulate \
  --m-input-file qiime2_results_sra/taxonomy.qza \
  --o-visualization qiime2_results_sra/taxonomy.qzv

qiime taxa barplot \
  --i-table qiime2_results_sra/table_final.qza \
  --i-taxonomy qiime2_results_sra/taxonomy.qza \
  --m-metadata-file metadata/GSE244391_metadata.tsv \
  --o-visualization qiime2_results_sra/taxa-bar-plots.qzv

echo "=== Complete: $(date) ==="