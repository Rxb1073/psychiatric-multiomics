#!/bin/bash
#SBATCH --job-name=trim_GSE102741
#SBATCH --time=70:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=logs/trim_%j.out
#SBATCH --error=logs/trim_%j.err

module purge
module load bear-apps/2023a
module load Trim_Galore/0.6.10-GCCcore-12.3.0

cd /rds/projects/e/elhamsak-dissertation-autism/RNA_seq/GSE102741 || exit 1
mkdir -p trimmed logs

echo "=== TrimGalore Started: $(date) ==="

# Paired-end processing
for R1 in fastq_files/*_1.fastq.gz; do
  SAMPLE=$(basename $R1 _1.fastq.gz)
  R2=fastq_files/${SAMPLE}_2.fastq.gz
  if [ -f "$R2" ]; then
    trim_galore \
      --quality 20 \
      --length 36 \
      --paired \
      --cores 4 \
      --output_dir trimmed \
      $R1 $R2
    echo "Done: $SAMPLE"
  fi
done

echo "=== TrimGalore Complete: $(date) ==="
echo "Trimmed files: $(ls trimmed/*_val_*.fq.gz 2>/dev/null | wc -l)"