#!/bin/bash
#SBATCH --job-name=dada2_769228_fix
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=logs/dada2_fix_%j.out
#SBATCH --error=logs/dada2_fix_%j.err

module purge
module load bluebear
module load bear-apps/2022b
module load QIIME2/2024.2

cd /rds/projects/e/elhamsak-dissertation-autism/Microbiome/PRJNA769228 || exit 1
export TMPDIR=/rds/projects/e/elhamsak-dissertation-autism/tmp
mkdir -p qiime2_results_fixed logs

# Find the demultiplexed paired-end artifact from the original run
DEMUX=qiime2_results/demux-paired-end.qza
echo "Using demux: $DEMUX"

echo "=== DADA2 paired-end (corrected trunc: f=280, r=230 -> overlap preserved): $(date) ==="
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs "$DEMUX" \
  --p-trunc-len-f 280 \
  --p-trunc-len-r 230 \
  --p-trim-left-f 17 \
  --p-trim-left-r 21 \
  --p-n-threads 16 \
  --o-table qiime2_results_fixed/table.qza \
  --o-representative-sequences qiime2_results_fixed/rep-seqs.qza \
  --o-denoising-stats qiime2_results_fixed/denoising-stats.qza \
  --verbose

qiime metadata tabulate \
  --m-input-file qiime2_results_fixed/denoising-stats.qza \
  --o-visualization qiime2_results_fixed/denoising-stats.qzv

echo "=== Complete: $(date) ==="