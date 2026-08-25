#!/bin/bash
#SBATCH --job-name=sexctrl
#SBATCH --time=4:00:00
#SBATCH --mem=180G
#SBATCH --cpus-per-task=4
#SBATCH --output=/rds/projects/e/elhamsak-dissertation-autism/Epigenomics/logs/sexctrl_%j.out
#SBATCH --error=/rds/projects/e/elhamsak-dissertation-autism/Epigenomics/logs/sexctrl_%j.err
set -euo pipefail
module purge
module load bear-apps/2023a/live
module load R/4.4.1-gfbf-2023a
module load R-bundle-Bioconductor/3.19-foss-2023a-R-4.4.1
Rscript /rds/projects/e/elhamsak-dissertation-autism/Epigenomics/epigenomics-ChIP-seq/GSE215991/sex_control.R
