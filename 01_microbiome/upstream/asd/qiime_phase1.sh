#!/bin/bash
#SBATCH --job-name=p1
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=/rds/projects/e/elhamsak-dissertation-autism/Microbiome/logs/%x_%j.out
#SBATCH --error=/rds/projects/e/elhamsak-dissertation-autism/Microbiome/logs/%x_%j.err
set -euo pipefail
DATASET=$1; RESULTS=$2; METADATA=$3
BASE=/rds/projects/e/elhamsak-dissertation-autism/Microbiome
cd "$BASE/$DATASET" || { echo "cannot cd to $BASE/$DATASET"; exit 1; }
module purge; module load bluebear; module load bear-apps/2022b; module load QIIME2/2024.2
export TMPDIR="$BASE/tmp/job_${SLURM_JOB_ID}"; mkdir -p "$TMPDIR" "$RESULTS/tree"
echo "===== $DATASET phase1 | results=$RESULTS | meta=$METADATA ====="

echo "----- [1/4] tree -----"
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences "$RESULTS/rep-seqs.qza" \
  --p-n-threads 16 \
  --o-alignment "$RESULTS/tree/aligned-rep-seqs.qza" \
  --o-masked-alignment "$RESULTS/tree/masked-aligned-rep-seqs.qza" \
  --o-tree "$RESULTS/tree/unrooted-tree.qza" \
  --o-rooted-tree "$RESULTS/tree/rooted-tree.qza"

echo "----- [2/4] table summarize -----"
qiime feature-table summarize \
  --i-table "$RESULTS/table.qza" \
  --m-sample-metadata-file "$METADATA" \
  --o-visualization "$RESULTS/table.qzv"

echo "----- [3/4] depth stats + id overlap -----"
qiime tools export --input-path "$RESULTS/table.qza" --output-path "$RESULTS/exported-table"
biom summarize-table -i "$RESULTS/exported-table/feature-table.biom" \
                     -o "$RESULTS/exported-table/biom-summary.txt"
sed -n '1,15p' "$RESULTS/exported-table/biom-summary.txt"
export BIOMFILE="$RESULTS/exported-table/feature-table.biom" METAFILE="$METADATA"
python3 - <<'PY' || echo "(id-check failed, non-fatal)"
import os, biom
t = biom.load_table(os.environ["BIOMFILE"])
tab = set(t.ids('sample')); meta=set()
with open(os.environ["METAFILE"]) as f:
    next(f, None)
    for line in f:
        if line.startswith('#'): continue
        s=line.split('\t')[0].strip()
        if s: meta.add(s)
print("table samples   :", len(tab))
print("metadata samples:", len(meta))
print("overlap         :", len(tab & meta))
print("table-only      :", len(tab - meta), sorted(tab-meta)[:10])
print("metadata-only   :", len(meta - tab), sorted(meta-tab)[:10])
PY

echo "----- [4/4] alpha-rarefaction -----"
MAXDEPTH=$(awk -F'[: ]+' '/Max:/{print int($3); exit}' "$RESULTS/exported-table/biom-summary.txt")
[ -n "${MAXDEPTH:-}" ] || { echo "could not parse max depth"; exit 1; }
echo "using --p-max-depth $MAXDEPTH"
qiime diversity alpha-rarefaction \
  --i-table "$RESULTS/table.qza" \
  --i-phylogeny "$RESULTS/tree/rooted-tree.qza" \
  --p-max-depth "$MAXDEPTH" \
  --m-metadata-file "$METADATA" \
  --o-visualization "$RESULTS/alpha-rarefaction.qzv"

echo "===== phase1 done: $DATASET ====="
ls -lh "$RESULTS"/tree/rooted-tree.qza "$RESULTS"/table.qzv "$RESULTS"/alpha-rarefaction.qzv
