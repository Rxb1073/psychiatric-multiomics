#!/bin/bash
#SBATCH --job-name=picrust2
#SBATCH --array=0-1
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --account=elhamsak-4-dissertat
#SBATCH --output=logs/picrust2_%A_%a.out
#SBATCH --error=logs/picrust2_%A_%a.err
module purge
module load bear-apps/2022b/live
module load Miniforge3/24.1.2-0
ENVDIR=/rds/projects/e/elhamsak-dissertation-neuropsychiatric2/MDD/microbiome/conda_envs/picrust2
eval "$(conda shell.bash hook)"; conda activate $ENVDIR

PICRUST=/rds/projects/e/elhamsak-dissertation-neuropsychiatric2/MDD/microbiome/picrust2
DATASETS=(PRJNA591924 PRJNA776170)
D=${DATASETS[$SLURM_ARRAY_TASK_ID]}
cd $PICRUST/$D

# remove any prior output (picrust2 refuses to overwrite)
rm -rf picrust2_output

picrust2_pipeline.py \
  -s exported/seqs/dna-sequences.fasta \
  -i exported/table/feature-table.biom \
  -o picrust2_output \
  -p 8 \
  --stratified \
  --verbose

echo "=== $D done. Outputs: ==="
ls -la picrust2_output/
