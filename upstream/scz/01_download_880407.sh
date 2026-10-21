#!/bin/bash
#SBATCH --job-name=dl-880407
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=08:00:00
#SBATCH --output=PRJNA880407/logs/download-%j.out

BASE=/rds/projects/e/elhamsak-dissertation-neuropsychiatric/Microbiome
DEST="$BASE/PRJNA880407"
RAW="$DEST/raw"
INFO="$DEST/metadata/ena_file_info.tsv"

mkdir -p "$RAW"
cd "$RAW"

tail -n +2 "$INFO" | while IFS=$'\t' read -r run samp title alias ftp md5 rest; do
  IFS=';' read -ra URLS <<< "$ftp"
  IFS=';' read -ra SUMS <<< "$md5"
  for i in "${!URLS[@]}"; do
    url="${URLS[$i]}"
    want="${SUMS[$i]}"
    file=$(basename "$url")
    if [ -f "$file" ]; then
      have=$(md5sum "$file" | cut -d' ' -f1)
      if [ "$have" == "$want" ]; then
        echo "[skip ] $file"
        continue
      fi
      echo "[redo ] $file (md5 mismatch)"
      rm -f "$file"
    fi
    for attempt in 1 2 3; do
      echo "[get  ] $file (attempt $attempt)"
      if wget -q -c "https://${url}" -O "$file"; then
        have=$(md5sum "$file" | cut -d' ' -f1)
        if [ "$have" == "$want" ]; then
          echo "[ok   ] $file"
          break
        fi
        echo "[bad  ] $file md5 mismatch"
        rm -f "$file"
      fi
      if [ "$attempt" -eq 3 ]; then echo "[FAIL ] $file"; fi
    done
  done
done

echo ""
echo "=========================================="
echo " files in raw/: $(ls -1 "$RAW"/*.fastq.gz 2>/dev/null | wc -l)  (expected 200)"
echo " total size   : $(du -sh "$RAW" | cut -f1)"
echo "=========================================="
