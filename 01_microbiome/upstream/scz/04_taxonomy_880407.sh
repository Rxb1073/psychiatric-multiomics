#!/bin/bash
#SBATCH --job-name=tax-880407
#SBATCH --cpus-per-task=4
#SBATCH --mem=120G
#SBATCH --time=12:00:00
#SBATCH --output=logs/taxonomy-%j.out
set -euo pipefail

module purge
module load bluebear
module load bear-apps/2022b
module load QIIME2/2024.2

BASE=/rds/projects/e/elhamsak-dissertation-neuropsychiatric/Microbiome/PRJNA880407
export TMPDIR=/rds/projects/e/elhamsak-dissertation-neuropsychiatric/tmp/job_${SLURM_JOB_ID}
mkdir -p "$TMPDIR"
cd "$BASE"

CLASSIFIER=/rds/projects/e/elhamsak-dissertation-autism/Microbiome/silva-138-99-515-806-pretrained.qza

echo ">> classifying against SILVA-138 V4 515/806 (OOM-safe) ..."
qiime feature-classifier classify-sklearn \
  --i-classifier "$CLASSIFIER" \
  --i-reads qiime2_v4/rep-seqs.qza \
  --p-n-jobs 1 \
  --p-reads-per-batch 1000 \
  --o-classification qiime2_v4/taxonomy.qza \
  --verbose

qiime metadata tabulate \
  --m-input-file qiime2_v4/taxonomy.qza \
  --o-visualization qiime2_v4/taxonomy.qzv

echo ">> taxonomy done."
