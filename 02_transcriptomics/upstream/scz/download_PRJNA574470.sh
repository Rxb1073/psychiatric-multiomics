#!/bin/bash
#SBATCH --job-name=PRJNA574470_dl
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --output=logs/dl_%j.out
#SBATCH --error=logs/dl_%j.err
module purge
module load bluebear
cd /rds/projects/e/elhamsak-dissertation-neuropsychiatric/RNA_seq/PRJNA574470 || exit 1
mkdir -p logs fastq_files metadata

echo "=== ENA Download (self-retrying): $(date) ==="
# (Re)fetch the file manifest from ENA
curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA574470&result=read_run&fields=run_accession,fastq_ftp,fastq_bytes,fastq_md5" \
  -o metadata/ena_file_info.tsv
cut -f2 metadata/ena_file_info.tsv | tail -n +2 | tr ';' '\n' | sed 's|^|ftp://|' > metadata/download_urls.txt
TOTAL=$(grep -c . metadata/download_urls.txt)
echo "Files to download: $TOTAL"

# Build checksum file
awk -F'\t' 'NR>1 {
  split($2, paths, ";"); split($4, hashes, ";");
  for(i in paths){ split(paths[i], p, "/"); print hashes[i]"  fastq_files/"p[length(p)] }
}' metadata/ena_file_info.tsv > metadata/checksums.md5

# Self-retrying loop: download, delete failures, repeat until clean (max 6 attempts)
for attempt in 1 2 3 4 5 6; do
  echo "--- Attempt $attempt ($(date)) ---"
  cd fastq_files
  cat ../metadata/download_urls.txt | xargs -n 1 -P 4 wget -q -c -T 60 -t 5
  cd ..
  # remove any file failing MD5 so it re-downloads next round
  md5sum -c metadata/checksums.md5 2>/dev/null | grep FAILED | sed 's/: FAILED.*//' | xargs -r rm -f
  GOOD=$(ls fastq_files/*.fastq.gz 2>/dev/null | wc -l)
  echo "  Good files so far: $GOOD / $TOTAL"
  if md5sum -c metadata/checksums.md5 --status 2>/dev/null; then
    echo "=== ALL FILES VERIFIED on attempt $attempt: $(date) ==="
    break
  fi
done
echo "=== Done: $(date) ==="
md5sum -c metadata/checksums.md5 2>/dev/null | grep -c FAILED | xargs echo "Remaining failures:"
