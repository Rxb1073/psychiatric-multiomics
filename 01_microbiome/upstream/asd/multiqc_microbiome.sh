#!/bin/bash
#SBATCH --job-name=multiqc_microbiome
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --output=/rds/projects/e/elhamsak-dissertation-autism/Microbiome/multiqc_%j.out
#SBATCH --error=/rds/projects/e/elhamsak-dissertation-autism/Microbiome/multiqc_%j.err

module purge
module load bear-apps/2024a
module load MultiQC/1.32-foss-2024a


cd /rds/projects/e/elhamsak-dissertation-autism/Microbiome || exit 1
mkdir -p MultiQC_reports

echo "=== MultiQC Started: $(date) ==="

for DATASET in PRJNA453621 PRJNA1022510 PRJNA529598 PRJNA769228 PRJNA484120; do
  echo "Processing $DATASET..."
  multiqc $DATASET/qc_results/fastqc/ \
    -o MultiQC_reports/$DATASET \
    -n ${DATASET}_multiqc_report \
    --force
done

echo "=== MultiQC Completed: $(date) ==="
echo "Reports in: MultiQC_reports/"