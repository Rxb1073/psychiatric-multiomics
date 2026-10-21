#!/bin/bash
#SBATCH --job-name=retrim_empty_PRJNA434002
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --array=0-5
#SBATCH --output=logs/retrim_empty_%A_%a.out
#SBATCH --error=logs/retrim_empty_%A_%a.err

module purge
module load bear-apps/2023a
module load Trim_Galore/0.6.10-GCCcore-12.3.0

cd /rds/projects/e/elhamsak-dissertation-autism/RNA_seq/PRJNA434002_bulk || exit 1
mkdir -p trimmed logs

# Build the list of empty samples fresh (size < 100 bytes), then split across 6 array tasks
mapfile -t EMPTY < <(
  while read s; do
    f="trimmed/${s}_1_val_1.fq.gz"
    if [ -f "$f" ] && [ $(stat -c%s "$f") -lt 100 ]; then echo "$s"; fi
  done < metadata/keep_samples.txt
)

TOTAL=${#EMPTY[@]}
echo "=== Task $SLURM_ARRAY_TASK_ID: $TOTAL empty samples total ==="

for ((i=SLURM_ARRAY_TASK_ID; i<TOTAL; i+=6)); do
  SAMPLE="${EMPTY[$i]}"
  R1=fastq_files/${SAMPLE}_1.fastq.gz
  R2=fastq_files/${SAMPLE}_2.fastq.gz

  # Delete the empty stubs so nothing is skipped
  rm -f trimmed/${SAMPLE}_1_val_1.fq.gz trimmed/${SAMPLE}_2_val_2.fq.gz

  if [ -f "$R1" ] && [ -f "$R2" ]; then
    echo "Trimming $SAMPLE ($(date))"
    trim_galore \
      --quality 20 \
      --length 36 \
      --paired \
      --cores 4 \
      --output_dir trimmed \
      $R1 $R2
    echo "Done: $SAMPLE ($(date))"
  else
    echo "MISSING RAW for $SAMPLE"
  fi
done

echo "=== Task $SLURM_ARRAY_TASK_ID complete: $(date) ==="