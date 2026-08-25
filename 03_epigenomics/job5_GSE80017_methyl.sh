#!/bin/bash
#SBATCH --job-name=GSE80017_Methyl
#SBATCH --account=elhamsak-dissertation-autism
#SBATCH --qos=bbdefault
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --output=/rds/projects/e/elhamsak-dissertation-autism/Epigenomics/logs/GSE80017_%j.out
#SBATCH --error=/rds/projects/e/elhamsak-dissertation-autism/Epigenomics/logs/GSE80017_%j.err

set -euo pipefail
echo "Job: ${SLURM_JOB_ID} | GSE80017 | Start: $(date)"
module purge
module load bear-apps/2023a/live
module load R/4.4.1-gfbf-2023a

Rscript /rds/projects/e/elhamsak-dissertation-autism/Epigenomics/epigenomics-methylation/GSE80017/run_analysis.R \
    2>&1 | tee /rds/projects/e/elhamsak-dissertation-autism/Epigenomics/epigenomics-methylation/GSE80017/logs/analysis.log

echo "GSE80017 COMPLETE: $(date)"
