#!/bin/bash
#SBATCH --job-name=fix_ids_PRJNA1022510
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=logs/fix_ids_%j.out
#SBATCH --error=logs/fix_ids_%j.err

module purge
module load bluebear

cd /rds/projects/e/elhamsak-dissertation-autism/Microbiome/PRJNA1022510/fastq_files || exit 1

echo "=== Fixing read IDs: $(date) ==="

# Fix sequences.fastq.gz - replace SRR26228689. with empty
zcat sequences.fastq.gz | \
  awk 'NR%4==1 {sub(/SRR[0-9]+\./, "", $1)} {print}' | \
  gzip > sequences_fixed.fastq.gz

# Fix barcodes.fastq.gz - replace SRR26228690. with empty  
zcat barcodes.fastq.gz | \
  awk 'NR%4==1 {sub(/SRR[0-9]+\./, "", $1)} {print}' | \
  gzip > barcodes_fixed.fastq.gz

# Replace original files
mv sequences_fixed.fastq.gz sequences.fastq.gz
mv barcodes_fixed.fastq.gz barcodes.fastq.gz

echo "=== Done: $(date) ==="

# Verify first few IDs match
echo "First sequence ID:"
zcat sequences.fastq.gz | head -1
echo "First barcode ID:"
zcat barcodes.fastq.gz | head -1