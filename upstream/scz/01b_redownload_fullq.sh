#!/bin/bash
#SBATCH --job-name=redl-880407
#SBATCH --cpus-per-task=6
#SBATCH --mem=16G
#SBATCH --time=10:00:00
#SBATCH --output=logs/redownload-%j.out
set -euo pipefail

module purge
module load bluebear
module load bear-apps/2022b
module load SRA-Toolkit/3.0.10-gompi-2022b

BASE=/rds/projects/e/elhamsak-dissertation-neuropsychiatric/Microbiome/PRJNA880407
RAW="$BASE/raw_fullq"
TMP="$BASE/sra_tmp"
INFO="$BASE/metadata/ena_file_info.tsv"
mkdir -p "$RAW" "$TMP"
cd "$TMP"

# list of run accessions from the ENA table (column 1)
RUNS=$(tail -n +2 "$INFO" | cut -f1)

for run in $RUNS; do
  # skip if both gzipped fastqs already exist
  if [ -f "$RAW/${run}_1.fastq.gz" ] && [ -f "$RAW/${run}_2.fastq.gz" ]; then
    echo "[skip ] $run"
    continue
  fi
  echo "[get  ] $run"
  prefetch "$run" --output-directory "$TMP" --max-size 50g >/dev/null 2>&1 || { echo "[FAIL prefetch] $run"; continue; }
  fasterq-dump "$run" --split-files --outdir "$TMP" -e 6 >/dev/null 2>&1 || { echo "[FAIL dump] $run"; continue; }
  # gzip and move both reads
  for r in 1 2; do
    if [ -f "$TMP/${run}_${r}.fastq" ]; then
      gzip -c "$TMP/${run}_${r}.fastq" > "$RAW/${run}_${r}.fastq.gz"
      rm -f "$TMP/${run}_${r}.fastq"
    fi
  done
  rm -rf "$TMP/$run"   # remove the .sra dir
  echo "[ok   ] $run"
done

echo ""
echo "=== full-quality FASTQs: $(ls -1 "$RAW"/*.fastq.gz 2>/dev/null | wc -l)  (expect 200) ==="
