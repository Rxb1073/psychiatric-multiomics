#!/bin/bash
#SBATCH --job-name=fastqc_PRJNA574470
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --output=logs/fastqc_%j.out
#SBATCH --error=logs/fastqc_%j.err
module purge
module load bluebear
module load FastQC/0.11.9-Java-11
cd /rds/projects/e/elhamsak-dissertation-neuropsychiatric/RNA_seq/PRJNA574470 || exit 1
mkdir -p qc_results/fastqc logs
echo "=== FastQC Started: $(date) ==="
TOTAL=$(find fastq_files -name "*.fastq.gz" | wc -l)
echo "Running FastQC on $TOTAL files..."
fastqc fastq_files/*.fastq.gz -o qc_results/fastqc -t 16
REPORTS=$(find qc_results/fastqc -name "*.html" | wc -l)
echo "=== FastQC Completed: $(date) — $REPORTS reports ==="
