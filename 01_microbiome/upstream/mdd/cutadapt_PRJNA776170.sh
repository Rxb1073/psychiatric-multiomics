#!/bin/bash
#SBATCH --job-name=cutadapt_PRJNA776170
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=logs/cutadapt_%j.out
#SBATCH --error=logs/cutadapt_%j.err
module purge
module load bluebear
module load bear-apps/2022b
module load QIIME2/2024.2
cd /rds/projects/e/elhamsak-dissertation-neuropsychiatric/Microbiome/PRJNA776170 || exit 1
export TMPDIR=/rds/projects/e/elhamsak-dissertation-neuropsychiatric/tmp
mkdir -p qiime2_v4 metadata logs $TMPDIR

# 1. Build paired-end manifest
echo "=== Building manifest: $(date) ==="
printf "sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n" > metadata/manifest.tsv
for r1 in fastq_files/*_1.fastq.gz; do
  sample=$(basename $r1 _1.fastq.gz)
  r2=fastq_files/${sample}_2.fastq.gz
  printf "%s\t%s/%s\t%s/%s\n" "$sample" "$PWD" "$r1" "$PWD" "$r2"
done >> metadata/manifest.tsv
echo "Manifest lines: $(wc -l < metadata/manifest.tsv)"

# 2. Import paired-end reads
echo "=== Import: $(date) ==="
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path metadata/manifest.tsv \
  --output-path qiime2_v4/demux-paired-end.qza \
  --input-format PairedEndFastqManifestPhred33V2

# 3. cutadapt V4 extraction (515F / 806R)
echo "=== cutadapt V4 extraction: $(date) ==="
qiime cutadapt trim-paired \
  --i-demultiplexed-sequences qiime2_v4/demux-paired-end.qza \
  --p-front-f GTGYCAGCMGCCGCGGTAA \
  --p-front-r GGACTACNVGGGTWTCTAAT \
  --p-discard-untrimmed \
  --p-cores 8 \
  --o-trimmed-sequences qiime2_v4/demux-v4.qza \
  --verbose

# 4. Summarize V4 reads
echo "=== Summarize V4: $(date) ==="
qiime demux summarize \
  --i-data qiime2_v4/demux-v4.qza \
  --o-visualization qiime2_v4/demux-v4.qzv

echo "=== Done: $(date) ==="
