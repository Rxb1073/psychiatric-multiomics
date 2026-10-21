#!/bin/bash
#SBATCH --job-name=fc_perBAM
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --array=0-51
#SBATCH --output=logs/fc_perbam_%A_%a.out
#SBATCH --error=logs/fc_perbam_%A_%a.err

module purge
module load bear-apps/2022a
module load Subread/2.0.4-GCC-11.3.0

GSE=/rds/projects/e/elhamsak-dissertation-autism/RNA_seq/GSE102741
GTF=/rds/projects/e/elhamsak-dissertation-autism/RNA_seq/reference/genome/Homo_sapiens.GRCh38.109.gtf
cd $GSE || exit 1
mkdir -p counts/perbam logs
TMP=/rds/projects/e/elhamsak-dissertation-autism/tmp/fc_one_${SLURM_ARRAY_TASK_ID}
mkdir -p $TMP

# Pick the ONE BAM for this task
mapfile -t ALL_BAMS < <(ls aligned/*/Aligned.sortedByCoord.out.bam | sort)
BAM="${ALL_BAMS[$SLURM_ARRAY_TASK_ID]}"
SAMPLE=$(basename $(dirname "$BAM"))

OUT="counts/perbam/${SAMPLE}.txt"
# Skip if already done (real output, not empty)
if [ -f "$OUT" ] && [ $(stat -c%s "$OUT") -gt 1000 ]; then
  echo "Skip (done): $SAMPLE"
  exit 0
fi

echo "=== featureCounts on $SAMPLE ($(date)) ==="
echo "BAM size: $(ls -lh $BAM | awk '{print $5}')"

featureCounts -T 4 -p --primary --countReadPairs -s 2 \
  --tmpDir $TMP \
  -a $GTF \
  -o "$OUT" \
  "$BAM"

echo "=== Done $SAMPLE ($(date)) ==="
rm -rf $TMP