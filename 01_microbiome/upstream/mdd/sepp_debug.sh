#!/bin/bash
#SBATCH --job-name=sepp_dbg
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --account=elhamsak-4-dissertat
#SBATCH --output=logs/sepp_dbg_%j.out
#SBATCH --error=logs/sepp_dbg_%j.err
set -e                      # STOP at first failure (no more misleading "ALL DONE")
module purge
module load bear-apps/2022b/live
module load QIIME2/2024.2
export NUMBA_CACHE_DIR=$TMPDIR/numba_cache; mkdir -p $NUMBA_CACHE_DIR

# give SEPP a big, writable scratch on the node (a common cause of failure is
# /tmp filling up or being too small for SEPP's intermediates)
export TMPDIR=${TMPDIR:-/tmp}
mkdir -p $TMPDIR/sepp_work
export SEPP_TEMP=$TMPDIR/sepp_work

PHY=/rds/projects/e/elhamsak-dissertation-neuropsychiatric2/MDD/microbiome/phylogenetic
REF=/rds/projects/e/elhamsak-dissertation-autism/Microbiome/reference_data/sepp-refs-silva-128.qza
cd $PHY
echo "TMPDIR=$TMPDIR"; df -h $TMPDIR

echo "START: $(date)"
qiime fragment-insertion sepp \
  --i-representative-sequences merged_rep_seqs.qza \
  --i-reference-database $REF \
  --p-threads 16 \
  --o-tree sepp_tree.qza \
  --o-placements sepp_placements.qza \
  --verbose
echo "SEPP DONE: $(date)"
