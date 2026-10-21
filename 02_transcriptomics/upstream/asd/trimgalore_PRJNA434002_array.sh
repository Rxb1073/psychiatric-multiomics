#!/bin/bash
#SBATCH --job-name=trim_PRJ_array
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --array=0-5
#SBATCH --output=logs/trim_array_%A_%a.out
#SBATCH --error=logs/trim_array_%A_%a.err

module purge
module load bear-apps/2023a
module load Trim_Galore/0.6.10-GCCcore-12.3.0

cd /rds/projects/e/elhamsak-dissertation-autism/RNA_seq/PRJNA434002_bulk || exit 1
mkdir -p trimmed logs

# Get all R1 files into an array, split into 6 groups by array task ID
mapfile -t ALL_R1 < <(ls fastq_files/*_1.fastq.gz | sort)
TOTAL=${#ALL_R1[@]}

echo "=== Array task $SLURM_ARRAY_TASK_ID started: $(date) ==="

for ((i=SLURM_ARRAY_TASK_ID; i<TOTAL; i+=6)); do
  R1="${ALL_R1[$i]}"
  SAMPLE=$(basename $R1 _1.fastq.gz)
  R2=fastq_files/${SAMPLE}_2.fastq.gz

  # Skip if already trimmed
  if [ -f "trimmed/${SAMPLE}_1_val_1.fq.gz" ]; then
    echo "Skipping $SAMPLE - already trimmed"
    continue
  fi

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

echo "=== Array task $SLURM_ARRAY_TASK_ID Complete: $(date) ==="