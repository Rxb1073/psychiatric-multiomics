#!/bin/bash
#SBATCH --job-name=install_mmuphin
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --account=elhamsak-4-dissertat
#SBATCH --output=logs/install_mmuphin_%j.out
#SBATCH --error=logs/install_mmuphin_%j.err
module purge
module load bear-apps/2022b/live
module load R/4.3.1-foss-2022b
export R_LIBS_USER=$HOME/R/mmuphin_lib
mkdir -p $R_LIBS_USER
cd /rds/projects/e/elhamsak-dissertation-neuropsychiatric2/MDD/microbiome/merged_genus
Rscript install_mmuphin.R
