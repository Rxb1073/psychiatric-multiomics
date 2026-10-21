#!/bin/bash
#SBATCH --job-name=fix_ancombc
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --account=elhamsak-4-dissertat
#SBATCH --output=logs/fix_ancombc_%j.out
#SBATCH --error=logs/fix_ancombc_%j.err
module purge
module load bear-apps/2022b/live
module load R/4.3.1-foss-2022b
export R_LIBS_USER=$HOME/R/mmuphin_lib
cd /rds/projects/e/elhamsak-dissertation-neuropsychiatric2/MDD/microbiome/merged_genus
Rscript fix_ancombc_deps.R
