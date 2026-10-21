#!/bin/bash
#SBATCH --job-name=sort_bams
#SBATCH --account=elhamsak-2-dissertat
#SBATCH --qos=bbdefault
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --array=1-1%8
#SBATCH --output=/rds/projects/e/elhamsak-dissertation-autism-2/MDD/transcriptomics/logs/sort_%A_%a.out
#SBATCH --error=/rds/projects/e/elhamsak-dissertation-autism-2/MDD/transcriptomics/logs/sort_%A_%a.err
module purge
module load bear-apps/2022a
module load SAMtools/1.16.1-GCC-11.3.0
A2=/rds/projects/e/elhamsak-dissertation-autism-2/MDD/transcriptomics
LIST=$A2/logs/unsorted_list.txt
BAM_IN=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $LIST)
[ -z "$BAM_IN" ] && { echo "no entry at $SLURM_ARRAY_TASK_ID"; exit 0; }
DIR=$(dirname "$BAM_IN")
OUT=$DIR/Aligned.sortedByCoord.out.bam
TMP=/rds/projects/e/elhamsak-dissertation-autism-2/MDD/tmp/sort_${SLURM_ARRAY_TASK_ID}_$$
mkdir -p $TMP
echo "sorting $BAM_IN"
samtools sort -@ 8 -m 4G -T $TMP/s -o "$OUT" "$BAM_IN"
if [ $(stat -c%s "$OUT") -gt 1000000 ]; then
  echo "OK: $(ls -lh $OUT | awk '{print $5}')"
  samtools quickcheck "$OUT" && rm -f "$BAM_IN" || { echo "quickcheck FAILED"; exit 1; }
else
  echo "SORT FAILED: $BAM_IN"; exit 1
fi
rm -rf $TMP
