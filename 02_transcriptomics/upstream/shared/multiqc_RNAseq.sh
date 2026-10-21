#!/bin/bash
#SBATCH --job-name=multiqc_RNAseq
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --output=/rds/projects/e/elhamsak-dissertation-autism/RNA_seq/logs/multiqc_%j.out
#SBATCH --error=/rds/projects/e/elhamsak-dissertation-autism/RNA_seq/logs/multiqc_%j.err

module purge
module load bear-apps/2024a
module load MultiQC/1.32-foss-2024a

cd /rds/projects/e/elhamsak-dissertation-autism/RNA_seq || exit 1
mkdir -p MultiQC_reports logs

echo "=== MultiQC Started: $(date) ==="

for DATASET in GSE102741 GSE211154 GSE59288_GSE51264 PRJNA434002_bulk; do
  echo "Processing $DATASET..."
  multiqc $DATASET/qc_results/fastqc/ \
    -o MultiQC_reports/$DATASET \
    -n ${DATASET}_multiqc_report \
    --force
done

echo "=== MultiQC Completed: $(date) ==="