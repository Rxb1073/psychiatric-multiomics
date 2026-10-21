#!/bin/bash
#SBATCH --job-name=fqc-880407
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=logs/fastqc-%j.out
set -euo pipefail

BASE=/rds/projects/e/elhamsak-dissertation-neuropsychiatric/Microbiome/PRJNA880407
cd "$BASE"
mkdir -p fastqc

echo ">> running FastQC on 200 files ..."
module purge
module load bluebear
module load bear-apps/2022b
module load FastQC/0.11.9-Java-11
fastqc -t 8 -o fastqc raw/*.fastq.gz

echo ">> aggregating with MultiQC ..."
module purge
module load bluebear
module load bear-apps/2024a
module load MultiQC/1.32-foss-2024a
cd fastqc
multiqc . -n multiqc_880407 -f

echo ">> done. Open fastqc/multiqc_880407.html"
