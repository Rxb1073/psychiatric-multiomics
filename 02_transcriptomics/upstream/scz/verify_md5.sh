#!/bin/bash
#SBATCH --job-name=md5_574470
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --output=logs/md5_%j.out
#SBATCH --error=logs/md5_%j.err
cd /rds/projects/e/elhamsak-dissertation-neuropsychiatric/RNA_seq/PRJNA574470 || exit 1
echo "=== md5 verify started: $(date) ==="
if md5sum -c metadata/checksums.md5 --quiet; then
  echo "ALL FILES VERIFIED"
else
  echo "FAILURES:"; md5sum -c metadata/checksums.md5 2>/dev/null | grep -v ': OK$'
fi
echo "=== done: $(date) ==="
