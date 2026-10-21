#!/bin/bash
#SBATCH --job-name=fc_heavy5
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --array=1-5
#SBATCH --output=logs/fcheavy_%A_%a.out
#SBATCH --error=logs/fcheavy_%A_%a.err
cd /rds/projects/e/elhamsak-dissertation-neuropsychiatric/RNA_seq/PRJNA574470 || exit 1
mkdir -p counts/perbam logs
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" metadata/heavy5.txt)
BAM=/rds/projects/e/elhamsak-dissertation-neuropsychiatric2/RNA_seq/PRJNA574470/aligned/${SAMPLE}/Aligned.sortedByCoord.out.bam
OUT=counts/perbam/${SAMPLE}.txt
if [ -f "$OUT" ] && [ $(stat -c%s "$OUT") -gt 1000 ]; then echo "Skip $SAMPLE"; exit 0; fi
WORK=/rds/projects/e/elhamsak-dissertation-neuropsychiatric2/RNA_seq/PRJNA574470/tmp/heavy_${SAMPLE}; rm -rf $WORK; mkdir -p $WORK
PRIM=$WORK/primary.bam

# step 1: filter secondary alignments (samtools)
module purge; module load bear-apps/2024a; module load SAMtools/1.21-GCC-13.3.0
echo "=== $SAMPLE filter: $(date) ==="
samtools view -@ 8 -b -F 0x100 "$BAM" -o "$PRIM"
echo "  filtered: $(ls -lh $PRIM | awk '{print $5}')"

# step 2: count (fresh module env to avoid toolchain conflict)
module purge; module load bear-apps/2022a; module load Subread/2.0.4-GCC-11.3.0
echo "=== $SAMPLE featureCounts: $(date) ==="
featureCounts -T 8 -p --countReadPairs -s 2 --primary --tmpDir $WORK -a /rds/projects/e/elhamsak-dissertation-autism/RNA_seq/reference/genome/Homo_sapiens.GRCh38.109.gtf -o $OUT "$PRIM"
rm -rf $WORK
echo "=== $SAMPLE done: $(date) ==="
