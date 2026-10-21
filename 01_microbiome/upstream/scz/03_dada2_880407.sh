#!/bin/bash
#SBATCH --job-name=dada2-880407
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=logs/dada2-%j.out
set -euo pipefail

module purge
module load bluebear
module load bear-apps/2022b
module load QIIME2/2024.2

BASE=/rds/projects/e/elhamsak-dissertation-neuropsychiatric/Microbiome/PRJNA880407
export TMPDIR=/rds/projects/e/elhamsak-dissertation-neuropsychiatric/tmp/job_${SLURM_JOB_ID}
mkdir -p "$TMPDIR"
cd "$BASE"
mkdir -p qiime2_v4

RAW="$BASE/raw_fullq"
MAN="$BASE/manifest.tsv"

# ---- build paired-end manifest from raw_fullq/ ----
echo -e "sample-id\tforward-absolute-filepath\treverse-absolute-filepath" > "$MAN"
for r1 in "$RAW"/*_1.fastq.gz; do
  run=$(basename "$r1" _1.fastq.gz)
  echo -e "${run}\t${RAW}/${run}_1.fastq.gz\t${RAW}/${run}_2.fastq.gz" >> "$MAN"
done
echo ">> manifest lines (incl header): $(wc -l < "$MAN")   (expect 101)"

# ---- import ----
if [ ! -f qiime2_v4/demux.qza ]; then
  echo ">> importing paired-end reads ..."
  qiime tools import \
    --type 'SampleData[PairedEndSequencesWithQuality]' \
    --input-path "$MAN" \
    --input-format PairedEndFastqManifestPhred33V2 \
    --output-path qiime2_v4/demux.qza
fi

# ---- DADA2 denoise (primers already removed; trunc 220/220; no trimming) ----
echo ">> DADA2 denoise-paired ..."
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs qiime2_v4/demux.qza \
  --p-trunc-len-f 200 \
  --p-trunc-len-r 200 \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-n-threads 16 \
  --o-table qiime2_v4/table.qza \
  --o-representative-sequences qiime2_v4/rep-seqs.qza \
  --o-denoising-stats qiime2_v4/denoising-stats.qza \
  --verbose

# ---- quick summaries ----
qiime metadata tabulate \
  --m-input-file qiime2_v4/denoising-stats.qza \
  --o-visualization qiime2_v4/denoising-stats.qzv
qiime feature-table summarize \
  --i-table qiime2_v4/table.qza \
  --m-sample-metadata-file metadata/metadata_qiime2.tsv \
  --o-visualization qiime2_v4/table.qzv

echo ">> DADA2 done. Check qiime2_v4/denoising-stats.qzv for retention."
