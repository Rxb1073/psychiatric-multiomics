#!/bin/bash
#SBATCH --job-name=div-880407
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=logs/diversity-%j.out
set -euo pipefail

module purge
module load bluebear
module load bear-apps/2022b
module load QIIME2/2024.2

BASE=/rds/projects/e/elhamsak-dissertation-neuropsychiatric/Microbiome/PRJNA880407
export TMPDIR=/rds/projects/e/elhamsak-dissertation-neuropsychiatric/tmp/job_${SLURM_JOB_ID}
mkdir -p "$TMPDIR"
cd "$BASE"
Q=qiime2_v4
META=metadata/metadata_qiime2.tsv
mkdir -p "$Q/tree"

# ---------------- SET AFTER PASS 1 ----------------
DEPTH=7243            # leave EMPTY on pass 1; set (e.g. 7200) on pass 2
MAXDEPTH=25000    # rarefaction x-axis cap (~median depth)
# --------------------------------------------------

# ---- Phase 1: phylogenetic tree ----
if [ ! -f "$Q/tree/rooted-tree.qza" ]; then
  echo ">> building tree ..."
  qiime phylogeny align-to-tree-mafft-fasttree \
    --i-sequences "$Q/rep-seqs.qza" \
    --p-n-threads 16 \
    --o-alignment "$Q/tree/aligned-rep-seqs.qza" \
    --o-masked-alignment "$Q/tree/masked-aligned-rep-seqs.qza" \
    --o-tree "$Q/tree/unrooted-tree.qza" \
    --o-rooted-tree "$Q/tree/rooted-tree.qza"
else
  echo ">> tree exists, skipping."
fi

# ---- table summary + rarefaction (to choose depth) ----
if [ ! -f "$Q/table.qzv" ]; then
  qiime feature-table summarize \
    --i-table "$Q/table.qza" \
    --m-sample-metadata-file "$META" \
    --o-visualization "$Q/table.qzv"
fi
if [ ! -f "$Q/alpha-rarefaction.qzv" ]; then
  qiime diversity alpha-rarefaction \
    --i-table "$Q/table.qza" \
    --i-phylogeny "$Q/tree/rooted-tree.qza" \
    --p-max-depth "$MAXDEPTH" \
    --m-metadata-file "$META" \
    --o-visualization "$Q/alpha-rarefaction.qzv"
fi

# ---- STOP if depth not chosen ----
if [ -z "$DEPTH" ]; then
  echo "============================================================"
  echo " PASS 1 done. Inspect table.qzv + alpha-rarefaction.qzv,"
  echo " set DEPTH (candidate ~7200 keeps all 100), re-submit."
  echo "============================================================"
  exit 0
fi

# ---- Phase 2: core metrics ----
CMR="$Q/core-metrics-results"
if [ ! -d "$CMR" ]; then
  echo ">> core-metrics at depth $DEPTH ..."
  qiime diversity core-metrics-phylogenetic \
    --i-phylogeny "$Q/tree/rooted-tree.qza" \
    --i-table "$Q/table.qza" \
    --p-sampling-depth "$DEPTH" \
    --m-metadata-file "$META" \
    --output-dir "$CMR"
fi

# ---- alpha significance ----
for A in shannon faith_pd observed_features evenness; do
  OUT="$CMR/${A}-significance.qzv"; [ -f "$OUT" ] && continue
  qiime diversity alpha-group-significance \
    --i-alpha-diversity "$CMR/${A}_vector.qza" \
    --m-metadata-file "$META" \
    --o-visualization "$OUT"
done

# ---- beta significance (PERMANOVA on diagnosis) ----
for B in unweighted_unifrac weighted_unifrac bray_curtis jaccard; do
  OUT="$CMR/${B}-diagnosis.qzv"; [ -f "$OUT" ] && continue
  qiime diversity beta-group-significance \
    --i-distance-matrix "$CMR/${B}_distance_matrix.qza" \
    --m-metadata-file "$META" \
    --m-metadata-column diagnosis \
    --p-pairwise \
    --o-visualization "$OUT"
done

echo ">> diversity done -> $CMR/"
