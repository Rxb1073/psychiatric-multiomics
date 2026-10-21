#!/bin/bash
#SBATCH --job-name=fc_PRJNA574470
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --array=1-78%20
#SBATCH --output=logs/fc_%A_%a.out
#SBATCH --error=logs/fc_%A_%a.err
module purge
module load bear-apps/2022a
module load Subread/2.0.4-GCC-11.3.0
cd /rds/projects/e/elhamsak-dissertation-neuropsychiatric/RNA_seq/PRJNA574470 || exit 1
mkdir -p counts/perbam logs
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" metadata/star_samples.txt)
SAMPLE=$(echo "$LINE" | cut -f1); LAYOUT=$(echo "$LINE" | cut -f2)
BAM=/rds/projects/e/elhamsak-dissertation-neuropsychiatric2/RNA_seq/PRJNA574470/aligned/${SAMPLE}/Aligned.sortedByCoord.out.bam
OUT=counts/perbam/${SAMPLE}.txt
if [ -f "$OUT" ] && [ $(stat -c%s "$OUT") -gt 1000 ]; then echo "Skip $SAMPLE"; exit 0; fi
TMP=/rds/projects/e/elhamsak-dissertation-neuropsychiatric2/RNA_seq/PRJNA574470/tmp/fc_${SAMPLE}_${SLURM_ARRAY_TASK_ID}; rm -rf $TMP; mkdir -p $TMP
echo "=== fc $SAMPLE ($LAYOUT): $(date) ==="
if [ "$LAYOUT" = "PE" ]; then
  featureCounts -T 8 -p --countReadPairs -s 2 --primary --tmpDir $TMP -a /rds/projects/e/elhamsak-dissertation-autism/RNA_seq/reference/genome/Homo_sapiens.GRCh38.109.gtf -o $OUT $BAM
else
  featureCounts -T 8 -s 2 --primary --tmpDir $TMP -a /rds/projects/e/elhamsak-dissertation-autism/RNA_seq/reference/genome/Homo_sapiens.GRCh38.109.gtf -o $OUT $BAM
fi
rm -rf $TMP
echo "=== $SAMPLE done: $(date) ==="
