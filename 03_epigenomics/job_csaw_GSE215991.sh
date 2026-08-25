#!/bin/bash
#SBATCH --job-name=GSE215991_csaw
#SBATCH --account=elhamsak-4-dissertat
#SBATCH --qos=bbdefault
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --output=/rds/projects/e/elhamsak-dissertation-neuropsychiatric2/Epigenomics/logs/csaw_%j.out
#SBATCH --error=/rds/projects/e/elhamsak-dissertation-neuropsychiatric2/Epigenomics/logs/csaw_%j.err

set -euo pipefail
echo "Job ${SLURM_JOB_ID} | GSE215991 csaw | start $(date)"
module purge
module load bear-apps/2023a/live
module load R-bundle-Bioconductor/3.19-foss-2023a-R-4.4.1
export TMPDIR="/rds/projects/e/elhamsak-dissertation-autism/Epigenomics/epigenomics-ChIP-seq/GSE215991/tmp_${SLURM_JOB_ID}"
mkdir -p "$TMPDIR"
Rscript /rds/projects/e/elhamsak-dissertation-neuropsychiatric2/Epigenomics/scripts/csaw_GSE215991.R
rm -rf "$TMPDIR"
echo "done $(date)"
