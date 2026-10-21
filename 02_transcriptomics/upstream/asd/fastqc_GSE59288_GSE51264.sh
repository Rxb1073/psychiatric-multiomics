#!/bin/bash
#SBATCH --job-name=fastqc_GSE59288_GSE51264
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --output=logs/fastqc_%j.out
#SBATCH --error=logs/fastqc_%j.err

module purge
module load bluebear
module load FastQC/0.11.9-Java-11

cd /rds/projects/e/elhamsak-dissertation-autism/RNA_seq/GSE59288_GSE51264 || exit 1
mkdir -p qc_results/fastqc

echo "=== FastQC Started: $(date) ==="
TOTAL_FILES=$(find fastq_files -name "*.fastq.gz" | wc -l)
echo "Running FastQC on $TOTAL_FILES files..."

fastqc fastq_files/*.fastq.gz \
  -o qc_results/fastqc \
  --extract \
  -t 8

REPORTS=$(find qc_results/fastqc -name "*.html" | wc -l)
echo "=== FastQC Completed: $(date) ==="
echo "Reports generated: $REPORTS"