#!/bin/bash
#SBATCH --job-name=maaslin2_mdd
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --account=elhamsak-4-dissertat
#SBATCH --output=logs/maaslin2_%j.out
#SBATCH --error=logs/maaslin2_%j.err
module purge
module load bear-apps/2022b/live
module load Miniforge3/24.1.2-0
ENVDIR=/rds/projects/e/elhamsak-dissertation-neuropsychiatric2/MDD/microbiome/conda_envs/ancombc
eval "$(conda shell.bash hook)"
conda activate $ENVDIR
cd /rds/projects/e/elhamsak-dissertation-neuropsychiatric2/MDD/microbiome/merged_genus
Rscript maaslin2_mdd.R
