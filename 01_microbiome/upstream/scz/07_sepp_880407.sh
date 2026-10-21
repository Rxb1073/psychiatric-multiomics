#!/bin/bash
#SBATCH --job-name=sepp-880407
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=80G
#SBATCH --output=logs/sepp-%j.out
#SBATCH --error=logs/sepp-%j.err
set -euo pipefail

module purge; module load bluebear; module load bear-apps/2022b; module load QIIME2/2024.2

BASE=/rds/projects/e/elhamsak-dissertation-neuropsychiatric/Microbiome/PRJNA880407
REF=/rds/projects/e/elhamsak-dissertation-neuropsychiatric/Microbiome/reference_data/sepp-refs-silva-128.qza
Q="$BASE/qiime2_v4"
export TMPDIR=/rds/projects/e/elhamsak-dissertation-neuropsychiatric/tmp/sepp_${SLURM_JOB_ID}
mkdir -p "$TMPDIR" "$Q/sepp" "$BASE/logs"
cd "$BASE"

echo "===== SCZ PRJNA880407 SEPP | $(date) ====="
echo "----- [1/3] fragment-insertion SEPP (slow) -----"
qiime fragment-insertion sepp \
  --i-representative-sequences "$Q/rep-seqs.qza" \
  --i-reference-database "$REF" \
  --p-threads 8 \
  --o-tree "$Q/sepp/insertion-tree.qza" \
  --o-placements "$Q/sepp/insertion-placements.qza"

echo "----- [2/3] filter-features to inserted fragments -----"
qiime fragment-insertion filter-features \
  --i-table "$Q/table.qza" \
  --i-tree "$Q/sepp/insertion-tree.qza" \
  --o-filtered-table "$Q/sepp/table-sepp-filtered.qza" \
  --o-removed-table "$Q/sepp/table-sepp-removed.qza"

echo "----- [3/3] midpoint-root -----"
qiime phylogeny midpoint-root \
  --i-tree "$Q/sepp/insertion-tree.qza" \
  --o-rooted-tree "$Q/sepp/rooted-tree-sepp.qza"

echo "----- retention: ASVs kept vs removed -----"
qiime tools export --input-path "$Q/sepp/table-sepp-filtered.qza" --output-path "$TMPDIR/kept" >/dev/null 2>&1
qiime tools export --input-path "$Q/sepp/table-sepp-removed.qza"  --output-path "$TMPDIR/rem"  >/dev/null 2>&1
biom summarize-table -i "$TMPDIR/kept/feature-table.biom" -o "$Q/sepp/kept-summary.txt"
biom summarize-table -i "$TMPDIR/rem/feature-table.biom"  -o "$Q/sepp/removed-summary.txt" 2>/dev/null || echo "(no removed features)"
echo "KEPT:";    grep -E "Num samples|Num observations|Total count" "$Q/sepp/kept-summary.txt"
echo "REMOVED:"; grep -E "Num observations|Total count" "$Q/sepp/removed-summary.txt" 2>/dev/null || echo "  none"
echo "===== SEPP done | $(date) ====="
