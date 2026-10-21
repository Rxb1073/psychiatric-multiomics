#!/bin/bash
#SBATCH --job-name=csaw_v2
#SBATCH --account=elhamsak-4-dissertat
#SBATCH --qos=bbdefault
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --output=/rds/projects/e/elhamsak-dissertation-neuropsychiatric2/Epigenomics/logs/csaw_v2_%j.out
#SBATCH --error=/rds/projects/e/elhamsak-dissertation-neuropsychiatric2/Epigenomics/logs/csaw_v2_%j.err
set -euo pipefail
module purge
module load bear-apps/2023a/live
module load R-bundle-Bioconductor/3.19-foss-2023a-R-4.4.1
Rscript /rds/projects/e/elhamsak-dissertation-neuropsychiatric2/Epigenomics/scripts/csaw_GSE215991_v2.R
