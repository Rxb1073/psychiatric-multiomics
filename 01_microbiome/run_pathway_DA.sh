#!/bin/bash
#SBATCH --job-name=pathway_DA
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --account=elhamsak-4-dissertat
#SBATCH --output=pathway_DA_%j.out
#SBATCH --error=pathway_DA_%j.err
module purge
module load bear-apps/2022b/live
module load Miniforge3/24.1.2-0
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /rds/projects/e/elhamsak-dissertation-neuropsychiatric2/MDD/microbiome/conda_envs/ancombc
cd /rds/projects/e/elhamsak-dissertation-neuropsychiatric2/MDD/microbiome/picrust2
Rscript picrust2_pathway_DA.R
