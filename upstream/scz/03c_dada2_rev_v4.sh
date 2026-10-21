#!/bin/bash
#SBATCH --job-name=dada2r-880407
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=logs/dada2rev-%j.out
set -euo pipefail

module purge
module load bluebear
module load bear-apps/2022b
module load QIIME2/2024.2
module load SeqKit/2.8.2

BASE=/rds/projects/e/elhamsak-dissertation-neuropsychiatric/Microbiome/PRJNA880407
export TMPDIR=/rds/projects/e/elhamsak-dissertation-neuropsychiatric/tmp/job_${SLURM_JOB_ID}
mkdir -p "$TMPDIR"
cd "$BASE"

RAW="$BASE/raw_fullq"
RC="$BASE/raw_rc"          # reverse-complemented R2 reads
mkdir -p "$RC"

# ---- reverse-complement every R2 into raw_rc/ as <run>.fastq.gz ----
echo ">> reverse-complementing R2 reads ..."
for r2 in "$RAW"/*_2.fastq.gz; do
  run=$(basename "$r2" _2.fastq.gz)
  out="$RC/${run}.fastq.gz"
  [ -f "$out" ] && continue
  # seqkit if available, else awk fallback
  if command -v seqkit >/dev/null 2>&1; then
    seqkit seq -r -p -t DNA "$r2" -o "$out"
  else
    zcat "$r2" | awk 'NR%4==1||NR%4==3{print} NR%4==2{cmd="rev <<<\""$0"\" | tr ACGTacgt TGCAtgca"; cmd|getline rc; close(cmd); print rc} NR%4==0{print}' | gzip > "$out"
  fi
done
echo ">> RC files: $(ls -1 "$RC"/*.fastq.gz | wc -l) (expect 100)"

# ---- single-end manifest on the RC (now V4, forward-oriented) reads ----
MAN="$BASE/manifest_rc.tsv"
echo -e "sample-id\tabsolute-filepath" > "$MAN"
for f in "$RC"/*.fastq.gz; do
  run=$(basename "$f" .fastq.gz)
  echo -e "${run}\t${f}" >> "$MAN"
done

if [ ! -f qiime2_v4/demux-rc.qza ]; then
  qiime tools import \
    --type 'SampleData[SequencesWithQuality]' \
    --input-path "$MAN" \
    --input-format SingleEndFastqManifestPhred33V2 \
    --output-path qiime2_v4/demux-rc.qza
fi

echo ">> DADA2 denoise-single on reverse-complemented V4 reads (trunc 200) ..."
qiime dada2 denoise-single \
  --i-demultiplexed-seqs qiime2_v4/demux-rc.qza \
  --p-trunc-len 200 \
  --p-trim-left 0 \
  --p-n-threads 16 \
  --o-table qiime2_v4/table.qza \
  --o-representative-sequences qiime2_v4/rep-seqs.qza \
  --o-denoising-stats qiime2_v4/denoising-stats.qza \
  --verbose

qiime metadata tabulate --m-input-file qiime2_v4/denoising-stats.qza --o-visualization qiime2_v4/denoising-stats.qzv
qiime feature-table summarize --i-table qiime2_v4/table.qza --m-sample-metadata-file metadata/metadata_qiime2.tsv --o-visualization qiime2_v4/table.qzv

echo ">> reverse-read V4 DADA2 done."
