#!/bin/bash
#SBATCH --job-name=trim_GSE102556
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --array=1-48%12
#SBATCH --output=logs/trim_%A_%a.out
#SBATCH --error=logs/trim_%A_%a.err
module purge; module load bear-apps/2023a; module load Trim_Galore/0.6.10-GCCcore-12.3.0
cd /rds/projects/e/elhamsak-dissertation-neuropsychiatric2/MDD/transcriptomics/GSE102556 || exit 1
mkdir -p trimmed logs metadata
ls fastq_files/*_1.fastq.gz | sed -E 's|.*/||; s/_1\.fastq\.gz//' > metadata/samples.txt
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" metadata/samples.txt)
[ -z "$SAMPLE" ] && { echo "no sample at index $SLURM_ARRAY_TASK_ID"; exit 1; }
if [ -f "trimmed/${SAMPLE}_1_val_1.fq.gz" ] && [ $(stat -c%s "trimmed/${SAMPLE}_1_val_1.fq.gz") -gt 1000 ]; then
  echo "skip $SAMPLE"; exit 0; fi
trim_galore --paired --quality 20 --length 36 --cores 4 -o trimmed \
  fastq_files/${SAMPLE}_1.fastq.gz fastq_files/${SAMPLE}_2.fastq.gz
