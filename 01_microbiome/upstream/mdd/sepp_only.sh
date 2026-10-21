#!/bin/bash
#SBATCH --job-name=sepp_mdd2
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=120G
#SBATCH --account=elhamsak-4-dissertat
#SBATCH --output=logs/sepp2_%j.out
#SBATCH --error=logs/sepp2_%j.err
module purge
module load bear-apps/2022b/live
module load QIIME2/2024.2
export NUMBA_CACHE_DIR=$TMPDIR/numba_cache; mkdir -p $NUMBA_CACHE_DIR

PHY=/rds/projects/e/elhamsak-dissertation-neuropsychiatric2/MDD/microbiome/phylogenetic
REF=/rds/projects/e/elhamsak-dissertation-autism/Microbiome/reference_data/sepp-refs-silva-128.qza
cd $PHY

echo "START: $(date)"
# steps 1-2 already done: merged_rep_seqs.qza + merged_asv_table.qza exist
qiime fragment-insertion sepp \
  --i-representative-sequences merged_rep_seqs.qza \
  --i-reference-database $REF \
  --p-threads 32 \
  --o-tree sepp_tree.qza \
  --o-placements sepp_placements.qza
echo "SEPP DONE: $(date)"

qiime fragment-insertion filter-features \
  --i-table merged_asv_table.qza \
  --i-tree sepp_tree.qza \
  --o-filtered-table merged_asv_inserted.qza \
  --o-removed-table merged_asv_removed.qza

qiime tools export --input-path merged_asv_inserted.qza --output-path tmp_biom
biom summarize-table -i tmp_biom/feature-table.biom > depth_summary.txt
head -30 depth_summary.txt
echo "ALL DONE: $(date)"
