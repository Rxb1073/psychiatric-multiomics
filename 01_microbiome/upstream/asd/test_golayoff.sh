#!/bin/bash
#SBATCH --job-name=test_golayoff
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=logs/test_golayoff_%j.out
#SBATCH --error=logs/test_golayoff_%j.err

module purge
module load bluebear
module load bear-apps/2022b
module load QIIME2/2024.2

cd /rds/projects/e/elhamsak-dissertation-autism/Microbiome/PRJNA1022510 || exit 1
export TMPDIR=/rds/projects/e/elhamsak-dissertation-autism/tmp

count_reads () {  # $1 = demux qza
  qiime demux summarize --i-data "$1" --o-visualization /tmp/sum.qzv >/dev/null 2>&1
  qiime tools export --input-path /tmp/sum.qzv --output-path /tmp/sumexp >/dev/null 2>&1
  awk -F',' 'NR>1{s+=$2} END{printf "TOTAL READS MAPPED: %d\n", s}' /tmp/sumexp/per-sample-fastq-counts.csv
  rm -rf /tmp/sum.qzv /tmp/sumexp
}

echo "===== GOLAY OFF, no rev-comp ====="
qiime demux emp-single \
  --i-seqs qiime2_results_sra/emp-single-end-sequences.qza \
  --m-barcodes-file metadata/GSE244391_metadata.tsv \
  --m-barcodes-column BarcodeSequence \
  --p-no-golay-error-correction \
  --o-per-sample-sequences qiime2_results_sra/demux_noG_fwd.qza \
  --o-error-correction-details /tmp/d1.qza --verbose 2>&1 | tail -2
count_reads qiime2_results_sra/demux_noG_fwd.qza

echo "===== GOLAY OFF, rev-comp-mapping ====="
qiime demux emp-single \
  --i-seqs qiime2_results_sra/emp-single-end-sequences.qza \
  --m-barcodes-file metadata/GSE244391_metadata.tsv \
  --m-barcodes-column BarcodeSequence \
  --p-no-golay-error-correction \
  --p-rev-comp-mapping-barcodes \
  --o-per-sample-sequences qiime2_results_sra/demux_noG_rcm.qza \
  --o-error-correction-details /tmp/d2.qza --verbose 2>&1 | tail -2
count_reads qiime2_results_sra/demux_noG_rcm.qza

echo "=== Done: $(date) ==="