#!/bin/bash
#SBATCH --job-name=trim_GSE211154
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=logs/trim_%j.out
#SBATCH --error=logs/trim_%j.err

module purge
module load bear-apps/2023a
module load Trim_Galore/0.6.10-GCCcore-12.3.0

cd /rds/projects/e/elhamsak-dissertation-autism/RNA_seq/GSE211154 || exit 1
mkdir -p trimmed logs

echo "=== TrimGalore Started: $(date) ==="

# Single-end, high adapter content
for R1 in fastq_files/*.fastq.gz; do
  SAMPLE=$(basename $R1 .fastq.gz)
  trim_galore \
    --quality 20 \
    --length 36 \
    --stringency 3 \
    --cores 4 \
    --output_dir trimmed \
    $R1
  echo "Done: $SAMPLE"
done

echo "=== TrimGalore Complete: $(date) ==="
echo "Trimmed files: $(ls trimmed/*_trimmed.fq.gz 2>/dev/null | wc -l)"