#!/bin/bash
#SBATCH --job-name=GSE74193_Methyl
#SBATCH --account=elhamsak-4-dissertat
#SBATCH --qos=bbdefault
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --output=/rds/projects/e/elhamsak-dissertation-neuropsychiatric2/Epigenomics/logs/GSE74193_%j.out
#SBATCH --error=/rds/projects/e/elhamsak-dissertation-neuropsychiatric2/Epigenomics/logs/GSE74193_%j.err

set -euo pipefail
echo "Job: ${SLURM_JOB_ID} | GSE74193 | Start: $(date)"
module purge
module load bear-apps/2023a/live
module load R/4.4.1-gfbf-2023a

Rscript /rds/projects/e/elhamsak-dissertation-neuropsychiatric2/Epigenomics/epigenomics-methylation/GSE74193/run_analysis.R \
    2>&1 | tee /rds/projects/e/elhamsak-dissertation-neuropsychiatric2/Epigenomics/epigenomics-methylation/GSE74193/logs/analysis.log

echo "GSE74193 COMPLETE: $(date)"
