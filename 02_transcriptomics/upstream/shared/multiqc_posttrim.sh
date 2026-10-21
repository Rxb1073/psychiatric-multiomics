#!/bin/bash
#SBATCH --job-name=mqc_posttrim
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --output=/rds/projects/e/elhamsak-dissertation-neuropsychiatric/RNA_seq/logs/mqc_posttrim_%j.out
#SBATCH --error=/rds/projects/e/elhamsak-dissertation-neuropsychiatric/RNA_seq/logs/mqc_posttrim_%j.err
module purge
module load bear-apps/2024a
module load MultiQC/1.32-foss-2024a
cd /rds/projects/e/elhamsak-dissertation-neuropsychiatric/RNA_seq || exit 1
mkdir -p MultiQC_reports logs
for DATASET in PRJNA319583 PRJNA314463 PRJNA574470; do
  multiqc $DATASET/qc_results/fastqc_posttrim/ -o MultiQC_reports/${DATASET}_posttrim -n ${DATASET}_posttrim_multiqc --force
done
echo "=== post-trim MultiQC done: $(date) ==="
