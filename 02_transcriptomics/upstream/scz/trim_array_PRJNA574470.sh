#!/bin/bash
#SBATCH --job-name=trim_PRJNA574470
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --array=1-79%20
#SBATCH --output=logs/trim_%A_%a.out
#SBATCH --error=logs/trim_%A_%a.err
module purge
module load bear-apps/2023a
module load Trim_Galore/0.6.10-GCCcore-12.3.0
cd /rds/projects/e/elhamsak-dissertation-neuropsychiatric/RNA_seq/PRJNA574470 || exit 1
mkdir -p trimmed logs

LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" metadata/trim_samples.txt)
SAMPLE=$(echo "$LINE" | cut -f1)
LAYOUT=$(echo "$LINE" | cut -f2)

if [ "$LAYOUT" = "PE" ]; then
  if [ -f "trimmed/${SAMPLE}_1_val_1.fq.gz" ] && [ $(stat -c%s "trimmed/${SAMPLE}_1_val_1.fq.gz") -gt 1000 ]; then
    echo "Skip PE $SAMPLE"; exit 0; fi
  echo "=== PE $SAMPLE start: $(date) ==="
  trim_galore --quality 20 --length 36 --paired --cores 4 --output_dir trimmed     fastq_files/${SAMPLE}_1.fastq.gz fastq_files/${SAMPLE}_2.fastq.gz
else
  if [ -f "trimmed/${SAMPLE}_trimmed.fq.gz" ] && [ $(stat -c%s "trimmed/${SAMPLE}_trimmed.fq.gz") -gt 1000 ]; then
    echo "Skip SE $SAMPLE"; exit 0; fi
  echo "=== SE $SAMPLE start: $(date) ==="
  trim_galore --quality 20 --length 36 --cores 4 --output_dir trimmed     fastq_files/${SAMPLE}.fastq.gz
fi
echo "=== $SAMPLE done: $(date) ==="
