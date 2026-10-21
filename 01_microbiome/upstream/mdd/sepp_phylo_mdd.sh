#!/bin/bash
#SBATCH --job-name=sepp_mdd
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --account=elhamsak-4-dissertat
#SBATCH --output=logs/sepp_mdd_%j.out
#SBATCH --error=logs/sepp_mdd_%j.err
module purge
module load bear-apps/2022b/live
module load QIIME2/2024.2
export NUMBA_CACHE_DIR=$TMPDIR/numba_cache; mkdir -p $NUMBA_CACHE_DIR

ROOT=/rds/projects/e/elhamsak-dissertation-neuropsychiatric2
MB=$ROOT/MDD/microbiome
PHY=$MB/phylogenetic
HARM=$ROOT/MDD/harmonization
REF=/rds/projects/e/elhamsak-dissertation-autism/Microbiome/reference_data/sepp-refs-silva-128.qza
cd $PHY
echo "REF=[$REF]"; ls -la $REF

# ---- 1. merge rep-seqs from both cohorts (ONE tree for all ASVs) ----
qiime feature-table merge-seqs \
  --i-data $MB/PRJNA591924/qiime2_v4/rep-seqs.qza \
  --i-data $MB/PRJNA776170/qiime2_v4/rep-seqs.qza \
  --o-merged-data merged_rep_seqs.qza

# ---- 2. merge ASV tables, filtered to harmonized samples ----
qiime feature-table filter-samples \
  --i-table $MB/PRJNA776170/qiime2_v4/table.qza \
  --m-metadata-file $HARM/PRJNA776170_harmonized.tsv \
  --o-filtered-table 776170_asv_68.qza
qiime feature-table filter-samples \
  --i-table $MB/PRJNA591924/qiime2_v4/table.qza \
  --m-metadata-file $HARM/PRJNA591924_harmonized.tsv \
  --o-filtered-table 591924_asv_90.qza
qiime feature-table merge \
  --i-tables 591924_asv_90.qza --i-tables 776170_asv_68.qza \
  --o-merged-table merged_asv_table.qza

# ---- 3. SEPP: all ASVs onto the SHARED SILVA-128 backbone ----
qiime fragment-insertion sepp \
  --i-representative-sequences merged_rep_seqs.qza \
  --i-reference-database $REF \
  --p-threads 16 \
  --o-tree sepp_tree.qza \
  --o-placements sepp_placements.qza

# ---- 4. drop fragments that failed insertion ----
qiime fragment-insertion filter-features \
  --i-table merged_asv_table.qza \
  --i-tree sepp_tree.qza \
  --o-filtered-table merged_asv_inserted.qza \
  --o-removed-table merged_asv_removed.qza

# ---- 5. read depths (for choosing rarefaction depth) ----
qiime tools export --input-path merged_asv_inserted.qza --output-path tmp_biom
biom summarize-table -i tmp_biom/feature-table.biom > depth_summary.txt
head -30 depth_summary.txt
echo "=== inserted: ==="; qiime tools export --input-path merged_asv_inserted.qza --output-path chk_ins 2>/dev/null
echo "=== removed (failed insertion): ==="; qiime tools export --input-path merged_asv_removed.qza --output-path chk_rem 2>/dev/null
