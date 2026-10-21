#!/bin/bash
#SBATCH --job-name=lodo
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=/rds/projects/e/elhamsak-dissertation-autism/Microbiome/logs/%x_%j.out
#SBATCH --error=/rds/projects/e/elhamsak-dissertation-autism/Microbiome/logs/%x_%j.err
module purge
module load bear-apps/2022b/live
module load Miniforge3/24.1.2-0
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /rds/projects/e/elhamsak-dissertation-autism/Microbiome/conda_envs/microbiome_r
cd /rds/projects/e/elhamsak-dissertation-autism/Microbiome/pooled
Rscript lodo_check.R
