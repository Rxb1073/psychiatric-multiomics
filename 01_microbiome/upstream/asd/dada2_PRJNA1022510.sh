#!/bin/bash
#SBATCH --job-name=dada2_PRJNA1022510
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=logs/dada2_%j.out
#SBATCH --error=logs/dada2_%j.err

module purge
module load bear-apps/2022b/live
module load QIIME2/2024.2

export TMPDIR=/rds/projects/e/elhamsak-dissertation-autism/tmp
export TEMP=$TMPDIR
export TMP=$TMPDIR
export R_TEMPDIR=$TMPDIR
mkdir -p $TMPDIR

cd /rds/projects/e/elhamsak-dissertation-autism/Microbiome/PRJNA1022510 || exit 1
mkdir -p qiime2_results logs

echo "=== QIIME2 EMP Import Started: $(date) ==="

qiime tools import \
  --type EMPSingleEndSequences \
  --input-path /rds/projects/e/elhamsak-dissertation-autism/Microbiome/PRJNA1022510/fastq_files \
  --output-path qiime2_results/emp-single-end.qza

echo "=== Demultiplexing: $(date) ==="

qiime demux emp-single \
  --i-seqs qiime2_results/emp-single-end.qza \
  --m-barcodes-file /rds/projects/e/elhamsak-dissertation-autism/Microbiome/PRJNA1022510/metadata/GSE244391_metadata.tsv \
  --m-barcodes-column BarcodeSequence \
  --p-no-golay-error-correction \
  --o-per-sample-sequences qiime2_results/demux.qza \
  --o-error-correction-details qiime2_results/demux-details.qza

qiime demux summarize \
  --i-data qiime2_results/demux.qza \
  --o-visualization qiime2_results/demux.qzv

echo "=== DADA2 Started: $(date) ==="

qiime dada2 denoise-single \
  --i-demultiplexed-seqs qiime2_results/demux.qza \
  --p-trunc-len 0 \
  --p-trim-left 19 \
  --p-max-ee 5 \
  --p-n-reads-learn 250000 \
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

CLASSIFIER=/rds/projects/e/elhamsak-dissertation-autism/Microbiome/silva-138-99-515-806-pretrained.qza


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
  --m-metadata-file /rds/projects/e/elhamsak-dissertation-autism/Microbiome/PRJNA1022510/metadata/GSE244391_metadata.tsv \
  --o-visualization qiime2_results/taxa-bar-plots.qzv

echo "=== Pipeline Complete: $(date) ==="