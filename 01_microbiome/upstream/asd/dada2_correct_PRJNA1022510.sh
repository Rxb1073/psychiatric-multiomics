#!/bin/bash
#SBATCH --job-name=dada2_correct_PRJNA1022510
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=logs/dada2_correct_%j.out
#SBATCH --error=logs/dada2_correct_%j.err

module purge
module load bluebear
module load bear-apps/2022b
module load QIIME2/2024.2

cd /rds/projects/e/elhamsak-dissertation-autism/Microbiome/PRJNA1022510 || exit 1
export TMPDIR=/rds/projects/e/elhamsak-dissertation-autism/tmp

echo "=== DADA2 on correct demux (golay off, fwd, 8.2M reads): $(date) ==="
qiime dada2 denoise-single \
  --i-demultiplexed-seqs qiime2_results_sra/demux_noG_fwd.qza \
  --p-trim-left 5 \
  --p-trunc-len 215 \
  --p-n-threads 16 \
  --o-table qiime2_results_sra/table_final.qza \
  --o-representative-sequences qiime2_results_sra/rep-seqs_final.qza \
  --o-denoising-stats qiime2_results_sra/denoising-stats_final.qza \
  --verbose

qiime metadata tabulate \
  --m-input-file qiime2_results_sra/denoising-stats_final.qza \
  --o-visualization qiime2_results_sra/denoising-stats_final.qzv

echo "=== Complete: $(date) ==="