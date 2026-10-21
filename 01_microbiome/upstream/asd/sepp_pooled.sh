#!/bin/bash
#SBATCH --job-name=sepp_pooled
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=240G
#SBATCH --output=/rds/projects/e/elhamsak-dissertation-autism/Microbiome/logs/%x_%j.out
#SBATCH --error=/rds/projects/e/elhamsak-dissertation-autism/Microbiome/logs/%x_%j.err
set -euo pipefail
BASE=/rds/projects/e/elhamsak-dissertation-autism/Microbiome
REF=$BASE/reference_data/sepp-refs-silva-128.qza

# keep temp on the project volume, NOT node-local /tmp (likely cause of previous failure)
export TMPDIR=$BASE/tmp/sepp_${SLURM_JOB_ID}
mkdir -p "$TMPDIR"
echo "TMPDIR=$TMPDIR"
df -h "$TMPDIR" | tail -1

cd "$BASE/pooled/unifrac"
module purge; module load bluebear; module load bear-apps/2022b; module load QIIME2/2024.2

qiime fragment-insertion sepp \
  --i-representative-sequences merged_rep_seqs.qza \
  --i-reference-database "$REF" \
  --p-threads 16 \
  --o-tree rooted-tree-sepp-pooled.qza \
  --o-placements placements-pooled.qza \
  --verbose

qiime fragment-insertion filter-features \
  --i-table merged_asv_table.qza \
  --i-tree rooted-tree-sepp-pooled.qza \
  --o-filtered-table merged_asv_inserted.qza \
  --o-removed-table merged_asv_removed.qza

echo "SEPP done $(date)"
rm -rf "$TMPDIR"
