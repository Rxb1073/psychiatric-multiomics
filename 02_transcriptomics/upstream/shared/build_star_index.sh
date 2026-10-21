#!/bin/bash
#SBATCH --job-name=star_index
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=/rds/projects/e/elhamsak-dissertation-autism/RNA_seq/logs/star_index_%j.out
#SBATCH --error=/rds/projects/e/elhamsak-dissertation-autism/RNA_seq/logs/star_index_%j.err

module purge
module load bear-apps/2024a
module load STAR/2.7.11b-GCC-13.3.0

mkdir -p /rds/projects/e/elhamsak-dissertation-autism/RNA_seq/reference/GRCh38
mkdir -p /rds/projects/e/elhamsak-dissertation-autism/RNA_seq/reference/genome

cd /rds/projects/e/elhamsak-dissertation-autism/RNA_seq/reference/genome

echo "=== Downloading GRCh38 reference: $(date) ==="

wget -c https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget -c https://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.gtf.gz

gunzip -k Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip -k Homo_sapiens.GRCh38.109.gtf.gz

echo "=== Building STAR index: $(date) ==="

STAR --runMode genomeGenerate \
  --genomeDir /rds/projects/e/elhamsak-dissertation-autism/RNA_seq/reference/GRCh38 \
  --genomeFastaFiles /rds/projects/e/elhamsak-dissertation-autism/RNA_seq/reference/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  --sjdbGTFfile /rds/projects/e/elhamsak-dissertation-autism/RNA_seq/reference/genome/Homo_sapiens.GRCh38.109.gtf \
  --runThreadN 16 \
  --sjdbOverhang 99

echo "=== STAR index complete: $(date) ==="
ls -lh /rds/projects/e/elhamsak-dissertation-autism/RNA_seq/reference/GRCh38/