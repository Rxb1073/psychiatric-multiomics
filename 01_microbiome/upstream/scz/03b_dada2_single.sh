#!/bin/bash
#SBATCH --job-name=dada2s-880407
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=logs/dada2single-%j.out
set -euo pipefail

module purge
module load bluebear
module load bear-apps/2022b
module load QIIME2/2024.2

BASE=/rds/projects/e/elhamsak-dissertation-neuropsychiatric/Microbiome/PRJNA880407
export TMPDIR=/rds/projects/e/elhamsak-dissertation-neuropsychiatric/tmp/job_${SLURM_JOB_ID}
mkdir -p "$TMPDIR"
cd "$BASE"

# import forward reads only (single-end) via a single-end manifest
RAW="$BASE/raw_fullq"
MAN="$BASE/manifest_single.tsv"
echo -e "sample-id\tabsolute-filepath" > "$MAN"
for r1 in "$RAW"/*_1.fastq.gz; do
  run=$(basename "$r1" _1.fastq.gz)
  echo -e "${run}\t${RAW}/${run}_1.fastq.gz" >> "$MAN"
done
echo ">> manifest lines: $(wc -l < "$MAN") (expect 101)"

if [ ! -f qiime2_v4/demux-single.qza ]; then
  qiime tools import \
    --type 'SampleData[SequencesWithQuality]' \
    --input-path "$MAN" \
    --input-format SingleEndFastqManifestPhred33V2 \
    --output-path qiime2_v4/demux-single.qza
fi

echo ">> DADA2 denoise-single (trunc 200, forward reads only) ..."
qiime dada2 denoise-single \
  --i-demultiplexed-seqs qiime2_v4/demux-single.qza \
  --p-trunc-len 200 \
  --p-trim-left 0 \
  --p-n-threads 16 \
  --o-table qiime2_v4/table.qza \
  --o-representative-sequences qiime2_v4/rep-seqs.qza \
  --o-denoising-stats qiime2_v4/denoising-stats.qza \
  --verbose

qiime metadata tabulate --m-input-file qiime2_v4/denoising-stats.qza --o-visualization qiime2_v4/denoising-stats.qzv
qiime feature-table summarize --i-table qiime2_v4/table.qza --m-sample-metadata-file metadata/metadata_qiime2.tsv --o-visualization qiime2_v4/table.qzv

echo ">> single-end DADA2 done."
