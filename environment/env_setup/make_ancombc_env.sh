#!/bin/bash
#SBATCH --job-name=make_ancombc_env
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --account=elhamsak-4-dissertat
#SBATCH --output=logs/make_env_%j.out
#SBATCH --error=logs/make_env_%j.err
module purge
module load bear-apps/2022b/live
module load Miniforge3/24.1.2-0
ENVDIR=/rds/projects/e/elhamsak-dissertation-neuropsychiatric2/MDD/microbiome/conda_envs/ancombc
export CONDA_PKGS_DIRS=/rds/projects/e/elhamsak-dissertation-neuropsychiatric2/MDD/microbiome/conda_envs/pkgs
mkdir -p $(dirname $ENVDIR) $CONDA_PKGS_DIRS
eval "$(conda shell.bash hook)"
conda create -y -p $ENVDIR -c conda-forge -c bioconda \
  jq bioconductor-ancombc bioconductor-treesummarizedexperiment r-base=4.3
conda activate $ENVDIR
Rscript -e 'suppressMessages({library(TreeSummarizedExperiment);library(ANCOMBC)}); cat("BOTH LOADED (conda): TRUE\n")'
