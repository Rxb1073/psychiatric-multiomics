#!/bin/bash
#SBATCH --job-name=seppdiv-880407
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/seppdiv-%j.out
#SBATCH --error=logs/seppdiv-%j.err
set -euo pipefail
module purge; module load bluebear; module load bear-apps/2022b; module load QIIME2/2024.2
BASE=/rds/projects/e/elhamsak-dissertation-neuropsychiatric/Microbiome/PRJNA880407
export TMPDIR=/rds/projects/e/elhamsak-dissertation-neuropsychiatric/tmp/seppdiv_${SLURM_JOB_ID}
mkdir -p "$TMPDIR"; cd "$BASE"
Q=qiime2_v4; META=metadata/metadata_qiime2.tsv
CM="$Q/core-metrics-sepp"
DEPTH=7243

rm -rf "$CM"
echo "----- core-metrics-phylogenetic on SEPP tree (depth $DEPTH) -----"
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny "$Q/sepp/rooted-tree-sepp.qza" \
  --i-table "$Q/sepp/table-sepp-filtered.qza" \
  --p-sampling-depth $DEPTH \
  --m-metadata-file "$META" \
  --output-dir "$CM"

echo "----- alpha significance -----"
for A in shannon faith_pd observed_features evenness; do
  qiime diversity alpha-group-significance \
    --i-alpha-diversity "$CM/${A}_vector.qza" \
    --m-metadata-file "$META" \
    --o-visualization "$CM/${A}-significance.qzv"
done

echo "----- beta significance: PERMANOVA + PERMDISP -----"
for B in unweighted_unifrac weighted_unifrac bray_curtis jaccard; do
  qiime diversity beta-group-significance \
    --i-distance-matrix "$CM/${B}_distance_matrix.qza" \
    --m-metadata-file "$META" --m-metadata-column diagnosis --p-pairwise \
    --o-visualization "$CM/${B}-diagnosis.qzv"
  qiime diversity beta-group-significance \
    --i-distance-matrix "$CM/${B}_distance_matrix.qza" \
    --m-metadata-file "$META" --m-metadata-column diagnosis --p-method permdisp \
    --o-visualization "$CM/${B}-permdisp.qzv"
done
echo ">> SEPP-tree diversity done -> $CM"
