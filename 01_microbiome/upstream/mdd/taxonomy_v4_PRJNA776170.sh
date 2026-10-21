#!/bin/bash
#SBATCH --job-name=tax_v4_PRJNA776170
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --output=logs/taxonomy_%j.out
#SBATCH --error=logs/taxonomy_%j.err
module purge
module load bluebear
module load bear-apps/2022b
module load QIIME2/2024.2
cd /rds/projects/e/elhamsak-dissertation-neuropsychiatric/Microbiome/PRJNA776170 || exit 1
export TMPDIR=/rds/projects/e/elhamsak-dissertation-neuropsychiatric/tmp
mkdir -p qiime2_v4 logs $TMPDIR
CLASSIFIER=/rds/projects/e/elhamsak-dissertation-autism/Microbiome/silva-138-99-515-806-pretrained.qza

echo "=== Taxonomy (OOM-safe: n-jobs 1, batch 1000): $(date) ==="
qiime feature-classifier classify-sklearn \
  --i-classifier $CLASSIFIER \
  --i-reads qiime2_v4/rep-seqs.qza \
  --p-n-jobs 1 \
  --p-reads-per-batch 1000 \
  --o-classification qiime2_v4/taxonomy.qza \
  --verbose

qiime metadata tabulate \
  --m-input-file qiime2_v4/taxonomy.qza \
  --o-visualization qiime2_v4/taxonomy.qzv

echo "=== Done: $(date) ==="