#!/bin/bash
#SBATCH --job-name=p2
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=/rds/projects/e/elhamsak-dissertation-autism/Microbiome/logs/%x_%j.out
#SBATCH --error=/rds/projects/e/elhamsak-dissertation-autism/Microbiome/logs/%x_%j.err
set -euo pipefail
DATASET=$1; RESULTS=$2; METADATA=$3; DEPTH=$4      # diagnosis col + Control ref are fixed below
BASE=/rds/projects/e/elhamsak-dissertation-autism/Microbiome
cd "$BASE/$DATASET" || { echo "cannot cd to $BASE/$DATASET"; exit 1; }
module purge; module load bluebear; module load bear-apps/2022b; module load QIIME2/2024.2
export TMPDIR="$BASE/tmp/job_${SLURM_JOB_ID}"; mkdir -p "$TMPDIR"
COL=diagnosis; REF=Control
echo "===== $DATASET phase2 | depth=$DEPTH | col=$COL ref=$REF ====="

CM="$RESULTS/core-metrics-results"
[ -d "$CM" ] && mv "$CM" "${CM}.bak.$(date +%s)"

echo "----- [1/4] core-metrics-phylogenetic -----"
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny "$RESULTS/tree/rooted-tree.qza" \
  --i-table "$RESULTS/table.qza" \
  --p-sampling-depth "$DEPTH" \
  --m-metadata-file "$METADATA" \
  --output-dir "$CM"

echo "----- [2/4] alpha-group-significance -----"
for m in shannon faith_pd observed_features evenness; do
  qiime diversity alpha-group-significance \
    --i-alpha-diversity "$CM/${m}_vector.qza" \
    --m-metadata-file "$METADATA" \
    --o-visualization "$CM/${m}-significance.qzv"
done

echo "----- [3/4] beta-group-significance (PERMANOVA) -----"
for dm in unweighted_unifrac weighted_unifrac bray_curtis jaccard; do
  qiime diversity beta-group-significance \
    --i-distance-matrix "$CM/${dm}_distance_matrix.qza" \
    --m-metadata-file "$METADATA" \
    --m-metadata-column "$COL" \
    --p-pairwise \
    --o-visualization "$CM/${dm}-${COL}.qzv"
done

echo "----- [4/4] ANCOM-BC at genus (un-rarefied table) -----"
qiime taxa collapse \
  --i-table "$RESULTS/table.qza" \
  --i-taxonomy "$RESULTS/taxonomy.qza" \
  --p-level 6 \
  --o-collapsed-table "$RESULTS/table-genus.qza"
qiime composition ancombc \
  --i-table "$RESULTS/table-genus.qza" \
  --m-metadata-file "$METADATA" \
  --p-formula "$COL" \
  --p-reference-levels "${COL}::${REF}" \
  --o-differentials "$RESULTS/ancombc-${COL}.qza"
qiime composition da-barplot \
  --i-data "$RESULTS/ancombc-${COL}.qza" \
  --p-significance-threshold 0.05 \
  --p-level-delimiter ';' \
  --o-visualization "$RESULTS/ancombc-${COL}.qzv"

echo "===== phase2 done: $DATASET ====="
ls -lh "$CM"/*-significance.qzv "$CM"/*-${COL}.qzv "$RESULTS/ancombc-${COL}.qzv"
