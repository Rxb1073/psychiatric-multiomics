#!/bin/bash
#SBATCH --job-name=fc_PRJNA319583
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --array=1-341%30
#SBATCH --output=logs/fc_%A_%a.out
#SBATCH --error=logs/fc_%A_%a.err
module purge
module load bear-apps/2022a
module load Subread/2.0.4-GCC-11.3.0
cd /rds/projects/e/elhamsak-dissertation-neuropsychiatric/RNA_seq/PRJNA319583 || exit 1
mkdir -p counts/perbam logs
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" metadata/star_samples.txt)
BAM=/rds/projects/e/elhamsak-dissertation-neuropsychiatric2/RNA_seq/PRJNA319583/aligned/${SAMPLE}/Aligned.sortedByCoord.out.bam
OUT=counts/perbam/${SAMPLE}.txt
# skip-logic: real output already there?
if [ -f "$OUT" ] && [ $(stat -c%s "$OUT") -gt 1000 ]; then echo "Skip $SAMPLE"; exit 0; fi
TMP=/rds/projects/e/elhamsak-dissertation-neuropsychiatric2/RNA_seq/PRJNA319583/tmp/fc_${SAMPLE}_${SLURM_ARRAY_TASK_ID}; rm -rf $TMP; mkdir -p $TMP
echo "=== fc $SAMPLE: $(date) ==="
featureCounts -T 4 -p --countReadPairs -s 0 --primary   --tmpDir $TMP   -a /rds/projects/e/elhamsak-dissertation-autism/RNA_seq/reference/genome/Homo_sapiens.GRCh38.109.gtf   -o $OUT   $BAM
rm -rf $TMP
echo "=== $SAMPLE done: $(date) ==="
