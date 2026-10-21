#!/bin/bash
#SBATCH --job-name=fastqc_GSE102556
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --output=logs/fastqc_%j.out
#SBATCH --error=logs/fastqc_%j.err
module purge; module load bluebear; module load FastQC/0.11.9-Java-11
cd /rds/projects/e/elhamsak-dissertation-neuropsychiatric2/MDD/transcriptomics/GSE102556 || exit 1
mkdir -p qc_results/fastqc
fastqc fastq_files/*.fastq.gz -o qc_results/fastqc --extract -t 8
echo "=== FastQC done: $(date) — $(ls qc_results/fastqc/*.html | wc -l) reports ==="
