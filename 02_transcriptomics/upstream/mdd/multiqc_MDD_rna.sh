#!/bin/bash
#SBATCH --job-name=multiqc_MDD_rna
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --output=/rds/projects/e/elhamsak-dissertation-neuropsychiatric2/MDD/transcriptomics/multiqc_%j.out
#SBATCH --error=/rds/projects/e/elhamsak-dissertation-neuropsychiatric2/MDD/transcriptomics/multiqc_%j.err
module purge; module load bear-apps/2024a; module load MultiQC/1.32-foss-2024a
cd /rds/projects/e/elhamsak-dissertation-neuropsychiatric2/MDD/transcriptomics || exit 1
for D in GSE102556 GSE101521; do
  multiqc $D/qc_results/fastqc/ -o $D/qc_results/multiqc -n ${D}_multiqc_report --force
done
echo "=== MultiQC done: $(date) ==="
