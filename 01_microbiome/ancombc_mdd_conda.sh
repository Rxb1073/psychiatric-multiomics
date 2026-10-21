#!/bin/bash
#SBATCH --job-name=ancombc_conda
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --account=elhamsak-4-dissertat
#SBATCH --output=logs/ancombc_conda_%j.out
#SBATCH --error=logs/ancombc_conda_%j.err
module purge
module load bear-apps/2022b/live
module load Miniforge3/24.1.2-0
ENVDIR=/rds/projects/e/elhamsak-dissertation-neuropsychiatric2/MDD/microbiome/conda_envs/ancombc
eval "$(conda shell.bash hook)"
conda activate $ENVDIR
cd /rds/projects/e/elhamsak-dissertation-neuropsychiatric2/MDD/microbiome/merged_genus
Rscript ancombc_mdd_conda.R
