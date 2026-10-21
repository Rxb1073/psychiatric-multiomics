#!/bin/bash
#SBATCH --job-name=test_orient
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=logs/test_orient_%j.out
#SBATCH --error=logs/test_orient_%j.err

module purge
module load bluebear
module load bear-apps/2022b
module load QIIME2/2024.2

cd /rds/projects/e/elhamsak-dissertation-autism/Microbiome/PRJNA1022510 || exit 1
export TMPDIR=/rds/projects/e/elhamsak-dissertation-autism/tmp

echo "===== OPTION A: rev-comp-mapping-barcodes ====="
qiime demux emp-single \
  --i-seqs qiime2_results_sra/emp-single-end-sequences.qza \
  --m-barcodes-file metadata/GSE244391_metadata.tsv \
  --m-barcodes-column BarcodeSequence \
  --p-golay-error-correction \
  --p-rev-comp-mapping-barcodes \
  --o-per-sample-sequences qiime2_results_sra/demux_A.qza \
  --o-error-correction-details qiime2_results_sra/demux_A_details.qza \
  --verbose 2>&1 | tail -8
echo "EXIT A: $?"

echo "===== OPTION B: rev-comp-barcodes ====="
qiime demux emp-single \
  --i-seqs qiime2_results_sra/emp-single-end-sequences.qza \
  --m-barcodes-file metadata/GSE244391_metadata.tsv \
  --m-barcodes-column BarcodeSequence \
  --p-golay-error-correction \
  --p-rev-comp-barcodes \
  --o-per-sample-sequences qiime2_results_sra/demux_B.qza \
  --o-error-correction-details qiime2_results_sra/demux_B_details.qza \
  --verbose 2>&1 | tail -8
echo "EXIT B: $?"

echo "===== Which produced output? ====="
ls -lh qiime2_results_sra/demux_A.qza qiime2_results_sra/demux_B.qza 2>/dev/null