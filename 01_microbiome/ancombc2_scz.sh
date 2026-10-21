#!/bin/bash
#SBATCH --job-name=ancombc2_scz
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=logs/ancombc2_scz_%j.out
#SBATCH --error=logs/ancombc2_scz_%j.err
module purge
module load bear-apps/2022b/live
module load Miniforge3/24.1.2-0
eval "$(conda shell.bash hook)"
conda activate /rds/projects/e/elhamsak-dissertation-autism/Microbiome/conda_envs/microbiome_r
cd /rds/projects/e/elhamsak-dissertation-neuropsychiatric/Microbiome/PRJNA880407
Rscript ancombc2_scz.R
