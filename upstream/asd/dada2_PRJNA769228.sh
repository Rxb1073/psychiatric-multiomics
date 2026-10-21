#!/bin/bash
#SBATCH --job-name=dada2_PRJNA769228
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --output=logs/dada2_%j.out
#SBATCH --error=logs/dada2_%j.err

module purge
module load bear-apps/2022b/live
module load QIIME2/2024.2

export TMPDIR=/rds/projects/e/elhamsak-dissertation-autism/tmp
mkdir -p $TMPDIR

cd /rds/projects/e/elhamsak-dissertation-autism/Microbiome/PRJNA769228 || exit 1
mkdir -p qiime2_results logs

CLASSIFIER=/rds/projects/e/elhamsak-dissertation-autism/Microbiome/silva-138-99-515-806-pretrained.qza

echo "=== QIIME2 Import Started: $(date) ==="

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path metadata/manifest.tsv \
  --output-path qiime2_results/demux-paired-end.qza \
  --input-format PairedEndFastqManifestPhred33V2

qiime demux summarize \
  --i-data qiime2_results/demux-paired-end.qza \
  --o-visualization qiime2_results/demux.qzv

echo "=== DADA2 Started: $(date) ==="

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs qiime2_results/demux-paired-end.qza \
  --p-trunc-len-f 230 \
  --p-trunc-len-r 200 \
  --p-trim-left-f 17 \
  --p-trim-left-r 21 \
  --p-n-threads 16 \
  --o-table qiime2_results/table.qza \
  --o-representative-sequences qiime2_results/rep-seqs.qza \
  --o-denoising-stats qiime2_results/denoising-stats.qza \
  --verbose

qiime feature-table summarize \
  --i-table qiime2_results/table.qza \
  --o-visualization qiime2_results/table.qzv

qiime feature-table tabulate-seqs \
  --i-data qiime2_results/rep-seqs.qza \
  --o-visualization qiime2_results/rep-seqs.qzv

qiime metadata tabulate \
  --m-input-file qiime2_results/denoising-stats.qza \
  --o-visualization qiime2_results/denoising-stats.qzv

echo "=== Taxonomy Classification: $(date) ==="

qiime feature-classifier classify-sklearn \
  --i-classifier $CLASSIFIER \
  --i-reads qiime2_results/rep-seqs.qza \
  --p-n-jobs 16 \
  --o-classification qiime2_results/taxonomy.qza

qiime metadata tabulate \
  --m-input-file qiime2_results/taxonomy.qza \
  --o-visualization qiime2_results/taxonomy.qzv

qiime taxa barplot \
  --i-table qiime2_results/table.qza \
  --i-taxonomy qiime2_results/taxonomy.qza \
  --m-metadata-file /rds/projects/e/elhamsak-dissertation-autism/Microbiome/PRJNA769228/metadata/metadata_qiime2.tsv \
  --o-visualization qiime2_results/taxa-bar-plots.qzv

echo "=== Pipeline Complete: $(date) ==="