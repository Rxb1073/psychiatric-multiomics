#!/bin/bash
#SBATCH --job-name=multiqc_PRJNA591924
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --output=logs/multiqc_%j.out
#SBATCH --error=logs/multiqc_%j.err
module purge
module load bear-apps/2024a
module load MultiQC/1.32-foss-2024a
cd /rds/projects/e/elhamsak-dissertation-neuropsychiatric2/MDD/microbiome/PRJNA591924 || exit 1
mkdir -p qc_results/multiqc
multiqc qc_results/fastqc/ -o qc_results/multiqc -n PRJNA591924_multiqc_report --force
echo "=== MultiQC done: $(date) ==="
