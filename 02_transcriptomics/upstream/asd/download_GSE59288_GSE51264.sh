#!/bin/bash
#SBATCH --job-name=GSE59288_dl
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=4G
#SBATCH --output=logs/dl_%j.out
#SBATCH --error=logs/dl_%j.err

module purge
module load bluebear

cd /rds/projects/e/elhamsak-dissertation-autism/RNA_seq/GSE59288_GSE51264 || exit 1
mkdir -p logs fastq_files metadata

echo "=== ENA Download Started: $(date) ==="

# GSE59288 - 34 ASD samples
curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA254971&result=read_run&fields=run_accession,fastq_ftp,fastq_bytes,fastq_md5" \
  -o metadata/ena_GSE59288.tsv

# GSE51264 - 38 control samples (human only)
curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA222268&result=read_run&fields=run_accession,fastq_ftp,fastq_bytes,fastq_md5" \
  -o metadata/ena_GSE51264.tsv

# Combine download URLs
cut -f2 metadata/ena_GSE59288.tsv | tail -n +2 | tr ';' '\n' | sed 's|^|ftp://|' > metadata/download_urls.txt
cut -f2 metadata/ena_GSE51264.tsv | tail -n +2 | tr ';' '\n' | sed 's|^|ftp://|' >> metadata/download_urls.txt

TOTAL_FILES=$(wc -l < metadata/download_urls.txt)
echo "Files to download: $TOTAL_FILES"

TOTAL_BYTES_1=$(cut -f3 metadata/ena_GSE59288.tsv | tail -n +2 | tr ';' '\n' | awk '{sum+=$1} END {print sum}')
TOTAL_BYTES_2=$(cut -f3 metadata/ena_GSE51264.tsv | tail -n +2 | tr ';' '\n' | awk '{sum+=$1} END {print sum}')
TOTAL_GB=$(echo "scale=2; ($TOTAL_BYTES_1+$TOTAL_BYTES_2)/1024/1024/1024" | bc)
echo "Expected size: ${TOTAL_GB} GB"

cd fastq_files
cat ../metadata/download_urls.txt | xargs -n 1 -P 8 wget -q -c -T 20 -t 3
cd ..

DOWNLOADED=$(ls -1 fastq_files/*.fastq.gz 2>/dev/null | wc -l)
echo "Downloaded: $DOWNLOADED / $TOTAL_FILES files"

if [ $DOWNLOADED -ne $TOTAL_FILES ]; then
  echo "WARNING: File count mismatch!"
  exit 1
fi

echo "=== Download Complete: $(date) ==="
echo "SUCCESS: All files downloaded and verified"