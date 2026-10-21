#!/bin/bash
#SBATCH --job-name=GSE102741_dl
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=4G
#SBATCH --output=logs/dl_%j.out
#SBATCH --error=logs/dl_%j.err

module purge
module load bluebear

cd /rds/projects/e/elhamsak-dissertation-autism/RNA_seq/GSE102741 || exit 1
mkdir -p logs fastq_files metadata

echo "=== ENA Download Started: $(date) ==="

curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA398545&result=read_run&fields=run_accession,fastq_ftp,fastq_bytes,fastq_md5" \
  -o metadata/ena_file_info.tsv

cut -f1 metadata/ena_file_info.tsv | tail -n +2 > metadata/SRR_Acc_List.txt
echo "Accessions: $(wc -l < metadata/SRR_Acc_List.txt)"

cut -f2 metadata/ena_file_info.tsv | tail -n +2 | tr ';' '\n' | sed 's|^|ftp://|' > metadata/download_urls.txt
TOTAL_FILES=$(wc -l < metadata/download_urls.txt)
echo "Files to download: $TOTAL_FILES"

TOTAL_BYTES=$(cut -f3 metadata/ena_file_info.tsv | tail -n +2 | tr ';' '\n' | awk '{sum+=$1} END {print sum}')
TOTAL_GB=$(echo "scale=2; $TOTAL_BYTES/1024/1024/1024" | bc)
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

echo "=== Verifying MD5 checksums ==="
awk -F'\t' 'NR>1 {
  split($2, paths, ";");
  split($4, hashes, ";");
  for(i in paths) {
    split(paths[i], parts, "/");
    filename=parts[length(parts)];
    print hashes[i] " fastq_files/" filename
  }
}' metadata/ena_file_info.tsv > metadata/checksums.md5

if md5sum -c metadata/checksums.md5 --status; then
  echo "ALL FILES VERIFIED"
else
  echo "WARNING: Some files corrupted!"
  md5sum -c metadata/checksums.md5 | grep "FAILED"
  exit 1
fi

echo "=== Download Complete: $(date) ==="
echo "SUCCESS: All files downloaded and verified"