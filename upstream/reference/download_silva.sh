#!/bin/bash
#SBATCH --job-name=silva_download
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --output=/rds/projects/e/elhamsak-dissertation-autism/Microbiome/silva_dl_%j.out
#SBATCH --error=/rds/projects/e/elhamsak-dissertation-autism/Microbiome/silva_dl_%j.err

echo "=== Downloading SILVA 138 V4 classifier: $(date) ==="

wget -c https://data.qiime2.org/2024.2/common/silva-138-99-515-806-nb-classifier.qza \
  -O /rds/projects/e/elhamsak-dissertation-autism/Microbiome/silva-138-99-515-806-nb-classifier.qza

SIZE=$(ls -lh /rds/projects/e/elhamsak-dissertation-autism/Microbiome/silva-138-99-seqs-515-806.qza | awk '{print $5}')
echo "Downloaded size: $SIZE"
echo "=== Done: $(date) ==="