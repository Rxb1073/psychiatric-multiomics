#!/bin/bash
#SBATCH --job-name=rebuild2_PRJNA1022510
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=logs/rebuild2_%j.out
#SBATCH --error=logs/rebuild2_%j.err

module purge
module load bluebear
module load bear-apps/2022b
module load QIIME2/2024.2

cd /rds/projects/e/elhamsak-dissertation-autism/Microbiome/PRJNA1022510 || exit 1
export TMPDIR=/rds/projects/e/elhamsak-dissertation-autism/tmp
mkdir -p $TMPDIR logs

# --- Stage fixed EMP files ---
# Replace barcodes with the ID-corrected version, gzip both into a clean import dir
rm -rf emp_import/emp_single_fixed
mkdir -p emp_import/emp_single_fixed
cp emp_import/emp_single/sequences.fastq      emp_import/emp_single_fixed/sequences.fastq
cp emp_import/emp_single/barcodes_fixed.fastq emp_import/emp_single_fixed/barcodes.fastq
gzip -f emp_import/emp_single_fixed/sequences.fastq
gzip -f emp_import/emp_single_fixed/barcodes.fastq

# Fresh output dir so nothing collides
rm -rf qiime2_results_sra
mkdir -p qiime2_results_sra

echo "=== Importing EMP single-end (ID-fixed, real-quality): $(date) ==="
qiime tools import \
  --type EMPSingleEndSequences \
  --input-path emp_import/emp_single_fixed \
  --output-path qiime2_results_sra/emp-single-end-sequences.qza

echo "=== Demultiplexing with golay correction: $(date) ==="
qiime demux emp-single \
  --i-seqs qiime2_results_sra/emp-single-end-sequences.qza \
  --m-barcodes-file metadata/GSE244391_metadata.tsv \
  --m-barcodes-column BarcodeSequence \
  --p-golay-error-correction \
  --o-per-sample-sequences qiime2_results_sra/demux.qza \
  --o-error-correction-details qiime2_results_sra/demux-details.qza \
  --verbose

qiime demux summarize \
  --i-data qiime2_results_sra/demux.qza \
  --o-visualization qiime2_results_sra/demux-summary.qzv

echo "=== DADA2 denoise (trim-left 5, trunc-len 215): $(date) ==="
qiime dada2 denoise-single \
  --i-demultiplexed-seqs qiime2_results_sra/demux.qza \
  --p-trim-left 5 \
  --p-trunc-len 215 \
  --p-n-threads 16 \
  --o-table qiime2_results_sra/table.qza \
  --o-representative-sequences qiime2_results_sra/rep-seqs.qza \
  --o-denoising-stats qiime2_results_sra/denoising-stats.qza \
  --verbose

qiime metadata tabulate \
  --m-input-file qiime2_results_sra/denoising-stats.qza \
  --o-visualization qiime2_results_sra/denoising-stats.qzv

echo "=== Complete: $(date) ==="