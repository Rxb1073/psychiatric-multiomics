#!/bin/bash
#SBATCH --job-name=dl_GSE41826
#SBATCH --account=elhamsak-4-dissertat
#SBATCH --qos=bbdefault
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=/rds/projects/e/elhamsak-dissertation-neuropsychiatric2/Epigenomics/logs/dl_GSE41826_%j.out
#SBATCH --error=/rds/projects/e/elhamsak-dissertation-neuropsychiatric2/Epigenomics/logs/dl_GSE41826_%j.err
set -euo pipefail

BASE=/rds/projects/e/elhamsak-dissertation-neuropsychiatric2/Epigenomics/epigenomics-methylation/GSE41826
echo "[BASE=$BASE]"
[ -n "$BASE" ] || { echo "BASE empty, aborting"; exit 1; }

cd "$BASE/idat" || { echo "idat dir missing - run folder setup first"; exit 1; }

wget -c "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE41826&format=file" -O GSE41826_RAW.tar
tar -xvf GSE41826_RAW.tar
gunzip -f *.idat.gz
rm -f GSE41826_RAW.tar

N=$(ls *.idat 2>/dev/null | wc -l)
echo "IDAT files present: $N   (expect EVEN, ~290 for the full ~145-array series)"
if [ "$N" -eq 0 ]; then echo "FAIL: no IDATs extracted"; exit 1; fi
if [ $(( N % 2 )) -ne 0 ]; then
  echo "WARNING: odd IDAT count -> an orphan Grn/Red exists. Remove it before minfi."
fi

echo "--- checking for orphan Grn/Red (should print nothing) ---"
for f in *_Grn.idat; do b=${f%_Grn.idat}; [ -e "${b}_Red.idat" ] || echo "ORPHAN (no Red): $f"; done
for f in *_Red.idat; do b=${f%_Red.idat}; [ -e "${b}_Grn.idat" ] || echo "ORPHAN (no Grn): $f"; done

echo "=== dl_GSE41826 done: $(date) ==="
