#!/bin/bash
#SBATCH --job-name=ancom-880407
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/ancombc-%j.out
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

# collapse to genus (level 6)
if [ ! -f "$Q/table-genus.qza" ]; then
  qiime taxa collapse \
    --i-table "$Q/table.qza" \
    --i-taxonomy "$Q/taxonomy.qza" \
    --p-level 6 \
    --o-collapsed-table "$Q/table-genus.qza"
fi

# ANCOM-BC, diagnosis as the formula (HC as reference by default alpha order)
qiime composition ancombc \
  --i-table "$Q/table-genus.qza" \
  --m-metadata-file "$META" \
  --p-formula 'diagnosis' \
  --o-differentials "$Q/ancombc-diagnosis.qza"

# barplot of significant genera (q < 0.05)
qiime composition da-barplot \
  --i-data "$Q/ancombc-diagnosis.qza" \
  --p-significance-threshold 0.05 \
  --o-visualization "$Q/ancombc-diagnosis.qzv"

# also tabulate the full differentials for reading
qiime tools export --input-path "$Q/ancombc-diagnosis.qza" --output-path "$Q/_ancombc"

echo ">> ANCOM-BC done. See $Q/ancombc-diagnosis.qzv and $Q/_ancombc/"
