#!/bin/bash
#SBATCH --job-name=dada2_v4_PRJNA776170
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=logs/dada2_v4_%j.out
#SBATCH --error=logs/dada2_v4_%j.err
module purge
module load bluebear
module load bear-apps/2022b
module load QIIME2/2024.2
cd /rds/projects/e/elhamsak-dissertation-neuropsychiatric/Microbiome/PRJNA776170 || exit 1
export TMPDIR=/rds/projects/e/elhamsak-dissertation-neuropsychiatric/tmp
mkdir -p qiime2_v4 logs

echo "=== DADA2 on V4-extracted reads: $(date) ==="
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs qiime2_v4/demux-v4.qza \
  --p-trunc-len-f 0 \
  --p-trunc-len-r 0 \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-n-threads 16 \
  --o-table qiime2_v4/table.qza \
  --o-representative-sequences qiime2_v4/rep-seqs.qza \
  --o-denoising-stats qiime2_v4/denoising-stats.qza \
  --verbose

qiime metadata tabulate \
  --m-input-file qiime2_v4/denoising-stats.qza \
  --o-visualization qiime2_v4/denoising-stats.qzv

qiime feature-table summarize \
  --i-table qiime2_v4/table.qza \
  --o-visualization qiime2_v4/table.qzv

echo "=== Done: $(date) ==="