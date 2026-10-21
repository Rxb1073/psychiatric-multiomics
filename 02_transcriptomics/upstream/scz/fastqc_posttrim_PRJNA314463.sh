#!/bin/bash
#SBATCH --job-name=fqpt_PRJNA314463
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --array=1-164%30
#SBATCH --output=logs/fqpt_%A_%a.out
#SBATCH --error=logs/fqpt_%A_%a.err
module purge
module load bluebear
module load FastQC/0.11.9-Java-11
cd /rds/projects/e/elhamsak-dissertation-neuropsychiatric/RNA_seq/PRJNA314463 || exit 1
mkdir -p qc_results/fastqc_posttrim logs
FQ=$(sed -n "${SLURM_ARRAY_TASK_ID}p" metadata/posttrim_fastqc_list.txt)
B=$(basename $FQ .fq.gz)
if [ -f "qc_results/fastqc_posttrim/${B}_fastqc.zip" ]; then echo "Skip $B"; exit 0; fi
echo "=== FastQC $B start: $(date) ==="
fastqc "$FQ" -o qc_results/fastqc_posttrim -t 2
echo "=== $B done: $(date) ==="
