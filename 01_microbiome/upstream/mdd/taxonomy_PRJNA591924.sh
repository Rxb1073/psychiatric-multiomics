#!/bin/bash
#SBATCH --job-name=tax_PRJNA591924
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --output=logs/taxonomy_%j.out
#SBATCH --error=logs/taxonomy_%j.err
set -euo pipefail
module purge
module load bear-apps/2022b/live
module load QIIME2/2024.2
cd /rds/projects/e/elhamsak-dissertation-neuropsychiatric2/MDD/microbiome/PRJNA591924 || exit 1
export TMPDIR=/rds/projects/e/elhamsak-dissertation-neuropsychiatric2/MDD/tmp
export NUMBA_CACHE_DIR=$TMPDIR/numba_cache; export MPLCONFIGDIR=$TMPDIR/mpl_cache
mkdir -p $NUMBA_CACHE_DIR $MPLCONFIGDIR
CLASSIFIER=/rds/projects/e/elhamsak-dissertation-autism/Microbiome/silva-138-99-515-806-pretrained.qza
qiime feature-classifier classify-sklearn \
  --i-classifier $CLASSIFIER --i-reads qiime2_v4/rep-seqs.qza \
  --p-n-jobs 1 --p-reads-per-batch 1000 \
  --o-classification qiime2_v4/taxonomy.qza --verbose
qiime metadata tabulate --m-input-file qiime2_v4/taxonomy.qza --o-visualization qiime2_v4/taxonomy.qzv
echo "=== taxonomy done: $(date) ==="
ls -lh qiime2_v4/taxonomy.qza
