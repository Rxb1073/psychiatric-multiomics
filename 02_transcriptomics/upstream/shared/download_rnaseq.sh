#!/bin/bash
#SBATCH --job-name=dl_rnaseq_np
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=logs/dl_rnaseq_%j.out
#SBATCH --error=logs/dl_rnaseq_%j.err

module purge
module load bear-apps/2022a
module load SRA-Toolkit/3.0.3-gompi-2022a

ROOT=/rds/projects/e/elhamsak-dissertation-neuropsychiatric
cd "$ROOT/RNA_seq" || exit 1
export TMPDIR=$ROOT/tmp
mkdir -p fastq_files/sra_tmp logs

for ACC in PRJNA319583 PRJNA314463 PRJNA574470; do
  LIST="metadata/${ACC}_SRR_list.txt"
  if [ ! -s "$LIST" ]; then
    echo "WARNING: $LIST missing/empty - skipping $ACC"; continue
  fi
  echo "=== $ACC : $(grep -c SRR "$LIST") runs ($(date)) ==="
  while read SRR; do
    SRR=$(echo "$SRR" | tr -d '\r' | xargs)
    case "$SRR" in SRR*) ;; *) continue ;; esac
    if ls fastq_files/${SRR}*.fastq.gz >/dev/null 2>&1; then
      echo "skip (exists): $SRR"; continue
    fi
    echo "Downloading $SRR ($(date))"
    prefetch "$SRR" -O fastq_files/sra_tmp
    fasterq-dump "fastq_files/sra_tmp/$SRR/$SRR.sra" --split-files -O fastq_files -e 8
    gzip fastq_files/${SRR}*.fastq
    rm -rf "fastq_files/sra_tmp/$SRR"
  done < "$LIST"
done
echo "=== ALL RNA-seq DOWNLOADS DONE: $(date) ==="
