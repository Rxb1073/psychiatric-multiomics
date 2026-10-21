#!/bin/bash
#SBATCH --job-name=sepp
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=80G
#SBATCH --output=/rds/projects/e/elhamsak-dissertation-autism/Microbiome/logs/%x_%j.out
#SBATCH --error=/rds/projects/e/elhamsak-dissertation-autism/Microbiome/logs/%x_%j.err
set -euo pipefail
DATASET=$1; RESULTS=$2
BASE=/rds/projects/e/elhamsak-dissertation-autism/Microbiome
REF="$BASE/reference_data/sepp-refs-silva-128.qza"
cd "$BASE/$DATASET" || { echo "cannot cd to $BASE/$DATASET"; exit 1; }
module purge; module load bluebear; module load bear-apps/2022b; module load QIIME2/2024.2
export TMPDIR="$BASE/tmp/sepp_${SLURM_JOB_ID}"; mkdir -p "$TMPDIR" "$RESULTS/sepp"
echo "===== $DATASET SEPP | results=$RESULTS | $(date) ====="
echo "----- [1/3] fragment-insertion SEPP (slow) -----"
qiime fragment-insertion sepp \
  --i-representative-sequences "$RESULTS/rep-seqs.qza" \
  --i-reference-database "$REF" \
  --p-threads 8 \
  --o-tree "$RESULTS/sepp/insertion-tree.qza" \
  --o-placements "$RESULTS/sepp/insertion-placements.qza"
echo "----- [2/3] filter-features to inserted fragments -----"
qiime fragment-insertion filter-features \
  --i-table "$RESULTS/table.qza" \
  --i-tree "$RESULTS/sepp/insertion-tree.qza" \
  --o-filtered-table "$RESULTS/sepp/table-sepp-filtered.qza" \
  --o-removed-table "$RESULTS/sepp/table-sepp-removed.qza"
echo "----- [3/3] midpoint-root -----"
qiime phylogeny midpoint-root \
  --i-tree "$RESULTS/sepp/insertion-tree.qza" \
  --o-rooted-tree "$RESULTS/sepp/rooted-tree-sepp.qza"
echo "----- retention: how many ASVs kept vs removed -----"
qiime tools export --input-path "$RESULTS/sepp/table-sepp-filtered.qza" --output-path "$TMPDIR/kept" >/dev/null 2>&1
qiime tools export --input-path "$RESULTS/sepp/table-sepp-removed.qza"  --output-path "$TMPDIR/rem"  >/dev/null 2>&1
biom summarize-table -i "$TMPDIR/kept/feature-table.biom" -o "$RESULTS/sepp/kept-summary.txt"
biom summarize-table -i "$TMPDIR/rem/feature-table.biom"  -o "$RESULTS/sepp/removed-summary.txt" 2>/dev/null || echo "(no removed features)"
echo "KEPT:";    grep -E "Num samples|Num observations|Total count" "$RESULTS/sepp/kept-summary.txt"
echo "REMOVED:"; grep -E "Num observations|Total count" "$RESULTS/sepp/removed-summary.txt" 2>/dev/null || echo "  none"
echo "===== SEPP done: $DATASET | $(date) ====="
ls -lh "$RESULTS/sepp/rooted-tree-sepp.qza" "$RESULTS/sepp/table-sepp-filtered.qza"
