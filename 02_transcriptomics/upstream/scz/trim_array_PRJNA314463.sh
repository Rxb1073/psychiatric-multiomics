#!/bin/bash
#SBATCH --job-name=trim_PRJNA314463
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --array=1-82%30
#SBATCH --output=logs/trim_%A_%a.out
#SBATCH --error=logs/trim_%A_%a.err
module purge
module load bear-apps/2023a
module load Trim_Galore/0.6.10-GCCcore-12.3.0
cd /rds/projects/e/elhamsak-dissertation-neuropsychiatric/RNA_seq/PRJNA314463 || exit 1
mkdir -p trimmed logs

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" metadata/trim_samples.txt)
R1=fastq_files/${SAMPLE}_1.fastq.gz
R2=fastq_files/${SAMPLE}_2.fastq.gz

# skip-logic (precaution #10): size>1000 = real output, not empty
if [ -f "trimmed/${SAMPLE}_1_val_1.fq.gz" ] && [ $(stat -c%s "trimmed/${SAMPLE}_1_val_1.fq.gz") -gt 1000 ]; then
  echo "Skip $SAMPLE - already trimmed"; exit 0
fi

echo "=== $SAMPLE start: $(date) ==="
trim_galore --quality 20 --length 36 --paired --cores 4 --output_dir trimmed $R1 $R2
echo "=== $SAMPLE done: $(date) ==="
