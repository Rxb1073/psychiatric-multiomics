#!/bin/bash
#SBATCH --job-name=fastqc_posttrim
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --output=/rds/projects/e/elhamsak-dissertation-autism/RNA_seq/MultiQC_reports/fastqc_posttrim_%j.out
#SBATCH --error=/rds/projects/e/elhamsak-dissertation-autism/RNA_seq/MultiQC_reports/fastqc_posttrim_%j.err

module purge
module load bluebear
module load bear-apps/2023a/live
module load FastQC/0.11.9-Java-11
module load bear-apps/2024a
module load MultiQC/1.32-foss-2024a

BASE=/rds/projects/e/elhamsak-dissertation-autism/RNA_seq

echo "=== Post-trim FastQC Started: $(date) ==="

# GSE102741 — paired-end
mkdir -p $BASE/GSE102741/QC/fastqc_posttrim
fastqc \
  $BASE/GSE102741/trimmed/*_val_1.fq.gz \
  $BASE/GSE102741/trimmed/*_val_2.fq.gz \
  --outdir $BASE/GSE102741/QC/fastqc_posttrim/ \
  --threads 16
echo "GSE102741 FastQC done: $(date)"

# GSE211154 — single-end
mkdir -p $BASE/GSE211154/QC/fastqc_posttrim
fastqc \
  $BASE/GSE211154/trimmed/*_trimmed.fq.gz \
  --outdir $BASE/GSE211154/QC/fastqc_posttrim/ \
  --threads 16
echo "GSE211154 FastQC done: $(date)"

# GSE59288_GSE51264 — single-end
mkdir -p $BASE/GSE59288_GSE51264/QC/fastqc_posttrim
fastqc \
  $BASE/GSE59288_GSE51264/trimmed/*_trimmed.fq.gz \
  --outdir $BASE/GSE59288_GSE51264/QC/fastqc_posttrim/ \
  --threads 16
echo "GSE59288_GSE51264 FastQC done: $(date)"

echo "=== All FastQC Complete: $(date) ==="
echo "=== Running MultiQC: $(date) ==="

# MultiQC for each dataset
multiqc $BASE/GSE102741/QC/fastqc_posttrim/ \
  -o $BASE/MultiQC_reports/GSE102741_fastqc_posttrim \
  -n GSE102741_fastqc_posttrim_multiqc \
  --force

multiqc $BASE/GSE211154/QC/fastqc_posttrim/ \
  -o $BASE/MultiQC_reports/GSE211154_fastqc_posttrim \
  -n GSE211154_fastqc_posttrim_multiqc \
  --force

multiqc $BASE/GSE59288_GSE51264/QC/fastqc_posttrim/ \
  -o $BASE/MultiQC_reports/GSE59288_GSE51264_fastqc_posttrim \
  -n GSE59288_GSE51264_fastqc_posttrim_multiqc \
  --force

echo "=== MultiQC Complete: $(date) ==="
echo "Reports saved to: $BASE/MultiQC_reports/"