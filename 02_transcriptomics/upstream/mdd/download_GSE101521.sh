#!/bin/bash
#SBATCH --job-name=GSE101521_dl
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=4G
#SBATCH --output=logs/dl_%j.out
#SBATCH --error=logs/dl_%j.err
module purge; module load bluebear
cd /rds/projects/e/elhamsak-dissertation-neuropsychiatric2/MDD/transcriptomics/GSE101521 || exit 1
mkdir -p logs fastq_files metadata
TOTAL=$(wc -l < metadata/download_urls.txt)
echo "=== Start $(date) — target $TOTAL files ==="
for attempt in 1 2 3 4 5 6 7 8; do
  cd fastq_files
  cat ../metadata/download_urls.txt | xargs -n 1 -P 8 wget -q -c -T 60 -t 5
  cd ..
  DL=$(ls fastq_files/*.fastq.gz 2>/dev/null | wc -l)
  echo "attempt $attempt: $DL / $TOTAL"
  [ "$DL" -ge "$TOTAL" ] && break
  sleep 10
done
awk -F'\t' 'NR>1{split($2,p,";");split($4,h,";");for(i in p){split(p[i],q,"/");print h[i]"  fastq_files/"q[length(q)]}}' \
  metadata/ena_file_info.tsv > metadata/checksums.md5
md5sum -c metadata/checksums.md5 --status && echo "ALL FILES VERIFIED" || { echo "SOME FAILED"; md5sum -c metadata/checksums.md5 | grep FAILED; }
echo "=== Done $(date): $(ls fastq_files/*.fastq.gz | wc -l)/$TOTAL ==="
