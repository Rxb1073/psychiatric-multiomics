#!/bin/bash
#SBATCH --job-name=cutadapt_769
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=logs/cutadapt_%j.out
#SBATCH --error=logs/cutadapt_%j.err
module purge
module load bluebear
module load bear-apps/2022b
module load QIIME2/2024.2
cd /rds/projects/e/elhamsak-dissertation-autism/Microbiome/PRJNA769228 || exit 1
export TMPDIR=/rds/projects/e/elhamsak-dissertation-autism/tmp
mkdir -p qiime2_v4 logs

echo "=== cutadapt V4 extraction from V3-V4 reads: $(date) ==="
qiime cutadapt trim-paired \
  --i-demultiplexed-sequences qiime2_results/demux-paired-end.qza \
  --p-front-f GTGYCAGCMGCCGCGGTAA \
  --p-front-r GGACTACNVGGGTWTCTAAT \
  --p-discard-untrimmed \
  --p-cores 8 \
  --o-trimmed-sequences qiime2_v4/demux-v4.qza \
  --verbose

echo "=== Summarize V4-trimmed reads: $(date) ==="
qiime demux summarize \
  --i-data qiime2_v4/demux-v4.qza \
  --o-visualization qiime2_v4/demux-v4.qzv

echo "=== Done: $(date) ==="
