#!/bin/bash
#SBATCH --job-name=sepp_fix
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --account=elhamsak-4-dissertat
#SBATCH --output=logs/sepp_fix_%j.out
#SBATCH --error=logs/sepp_fix_%j.err
set -e
module purge
module load bear-apps/2022b/live
module load QIIME2/2024.2

PHY=/rds/projects/e/elhamsak-dissertation-neuropsychiatric2/MDD/microbiome/phylogenetic
REF=/rds/projects/e/elhamsak-dissertation-autism/Microbiome/reference_data/sepp-refs-silva-128.qza

# CRITICAL: TMPDIR on the project volume (terabytes), not the 4GB node /tmp.
# QIIME2 stages artifact data into TMPDIR *and* SEPP writes intermediates there.
export TMPDIR=$PHY/scratch_$SLURM_JOB_ID
mkdir -p $TMPDIR
export SEPP_TEMP=$TMPDIR
export NUMBA_CACHE_DIR=$TMPDIR/numba_cache; mkdir -p $NUMBA_CACHE_DIR
echo "TMPDIR=$TMPDIR"
df -h $TMPDIR

cd $PHY
echo "START: $(date)"
qiime fragment-insertion sepp \
  --i-representative-sequences merged_rep_seqs.qza \
  --i-reference-database $REF \
  --p-threads 16 \
  --o-tree sepp_tree.qza \
  --o-placements sepp_placements.qza \
  --verbose
echo "SEPP DONE: $(date)"
echo "scratch used:"; du -sh $TMPDIR
