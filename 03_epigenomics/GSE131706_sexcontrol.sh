#!/bin/bash
#SBATCH --time=00:20:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --output=sexcontrol_ASD.out
module purge
module load bluebear
module load bear-apps/2023a/live
module load R/4.4.1-gfbf-2023a
module load R-bundle-Bioconductor/3.19-foss-2023a-R-4.4.1
Rscript sexcontrol_ASD.R
