#!/bin/bash
#SBATCH --job-name=dada2_PRJNA591924
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=logs/dada2_%j.out
#SBATCH --error=logs/dada2_%j.err
module purge
module load bluebear
module load bear-apps/2022b
module load QIIME2/2024.2
cd /rds/projects/e/elhamsak-dissertation-neuropsychiatric2/MDD/microbiome/PRJNA591924 || exit 1
export TMPDIR=/rds/projects/e/elhamsak-dissertation-neuropsychiatric2/MDD/tmp
mkdir -p qiime2_v4 metadata logs $TMPDIR

echo "=== Manifest: $(date) ==="
printf "sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n" > metadata/manifest.tsv
for r1 in fastq_files/*_1.fastq.gz; do
  s=$(basename $r1 _1.fastq.gz); r2=fastq_files/${s}_2.fastq.gz
  printf "%s\t%s/%s\t%s/%s\n" "$s" "$PWD" "$r1" "$PWD" "$r2"
done >> metadata/manifest.tsv
echo "samples in manifest: $(( $(wc -l < metadata/manifest.tsv) - 1 ))"

echo "=== Import: $(date) ==="
qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path metadata/manifest.tsv \
  --output-path qiime2_v4/demux-paired-end.qza \
  --input-format PairedEndFastqManifestPhred33V2
qiime demux summarize --i-data qiime2_v4/demux-paired-end.qza \
  --o-visualization qiime2_v4/demux.qzv

echo "=== DADA2 (V4, primers pre-removed): $(date) ==="
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs qiime2_v4/demux-paired-end.qza \
  --p-trunc-len-f 230 --p-trunc-len-r 180 \
  --p-trim-left-f 0 --p-trim-left-r 0 \
  --p-n-threads 16 \
  --o-table qiime2_v4/table.qza \
  --o-representative-sequences qiime2_v4/rep-seqs.qza \
  --o-denoising-stats qiime2_v4/denoising-stats.qza --verbose

qiime metadata tabulate --m-input-file qiime2_v4/denoising-stats.qza \
  --o-visualization qiime2_v4/denoising-stats.qzv
qiime feature-table summarize --i-table qiime2_v4/table.qza \
  --o-visualization qiime2_v4/table.qzv
qiime feature-table tabulate-seqs --i-data qiime2_v4/rep-seqs.qza \
  --o-visualization qiime2_v4/rep-seqs.qzv
echo "=== DADA2 done: $(date) ==="
