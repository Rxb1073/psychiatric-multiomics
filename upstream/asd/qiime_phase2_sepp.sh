#!/bin/bash
#SBATCH --job-name=p2sepp
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=/rds/projects/e/elhamsak-dissertation-autism/Microbiome/logs/%x_%j.out
#SBATCH --error=/rds/projects/e/elhamsak-dissertation-autism/Microbiome/logs/%x_%j.err
set -euo pipefail
DATASET=$1; RESULTS=$2; METADATA=$3; DEPTH=$4
BASE=/rds/projects/e/elhamsak-dissertation-autism/Microbiome
cd "$BASE/$DATASET" || { echo "cannot cd to $BASE/$DATASET"; exit 1; }
module purge; module load bluebear; module load bear-apps/2022b; module load QIIME2/2024.2
export TMPDIR="$BASE/tmp/p2sepp_${SLURM_JOB_ID}"; mkdir -p "$TMPDIR"
COL=diagnosis
echo "===== $DATASET phase2-SEPP | depth=$DEPTH | $(date) ====="

CM="$RESULTS/core-metrics-sepp"
[ -d "$CM" ] && mv "$CM" "${CM}.bak.$(date +%s)"

echo "----- core-metrics-phylogenetic on SEPP tree -----"
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny "$RESULTS/sepp/rooted-tree-sepp.qza" \
  --i-table "$RESULTS/table.qza" \
  --p-sampling-depth "$DEPTH" \
  --m-metadata-file "$METADATA" \
  --output-dir "$CM"

echo "----- alpha significance (phylo: faith_pd; plus others for completeness) -----"
for m in faith_pd shannon observed_features evenness; do
  qiime diversity alpha-group-significance \
    --i-alpha-diversity "$CM/${m}_vector.qza" \
    --m-metadata-file "$METADATA" \
    --o-visualization "$CM/${m}-significance.qzv"
done

echo "----- beta significance (phylo UniFrac; plus bray/jaccard) -----"
for dm in unweighted_unifrac weighted_unifrac bray_curtis jaccard; do
  qiime diversity beta-group-significance \
    --i-distance-matrix "$CM/${dm}_distance_matrix.qza" \
    --m-metadata-file "$METADATA" \
    --m-metadata-column "$COL" \
    --p-pairwise \
    --o-visualization "$CM/${dm}-${COL}.qzv"
done

echo "===== phase2-SEPP done: $DATASET | $(date) ====="
ls -lh "$CM"/*-significance.qzv "$CM"/*-${COL}.qzv
