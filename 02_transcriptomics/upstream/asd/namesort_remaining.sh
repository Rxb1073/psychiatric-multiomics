#!/bin/bash
#SBATCH --job-name=namesort_rem
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --array=0-5
#SBATCH --output=logs/namesort_rem_%A_%a.out
#SBATCH --error=logs/namesort_rem_%A_%a.err
module purge
module load bear-apps/2022a
module load SAMtools/1.16.1-GCC-11.3.0
module load Subread/2.0.4-GCC-11.3.0
GSE=/rds/projects/e/elhamsak-dissertation-autism/RNA_seq/GSE102741
GTF=/rds/projects/e/elhamsak-dissertation-autism/RNA_seq/reference/genome/Homo_sapiens.GRCh38.109.gtf
cd $GSE || exit 1
mkdir -p counts/perbam namesorted logs
TMP=/rds/projects/e/elhamsak-dissertation-autism/tmp/ns_${SLURM_ARRAY_TASK_ID}
mkdir -p $TMP

# List of the 6 remaining samples (FILL IN from the check above)
SAMPLES=(SRR5938444 SRR5938445 SRR5938447 SRR5938466 SRR5938468 SRR5938469)
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

OUT=counts/perbam/${SAMPLE}.txt
if [ -f "$OUT" ] && [ $(stat -c%s "$OUT") -gt 1000 ]; then echo "Skip $SAMPLE"; exit 0; fi

BAM=$GSE/aligned/$SAMPLE/Aligned.sortedByCoord.out.bam
echo "=== Name-sort $SAMPLE ($(date)) size $(ls -lh $BAM|awk '{print $5}') ==="
samtools sort -n -@ 8 -m 4G -T $TMP/ns_$SAMPLE -o namesorted/${SAMPLE}.ns.bam "$BAM"
echo "=== featureCounts $SAMPLE ($(date)) ==="
featureCounts -T 8 -p --countReadPairs -s 2 --tmpDir $TMP \
  -a $GTF -o "$OUT" namesorted/${SAMPLE}.ns.bam
echo "=== Done $SAMPLE ($(date)) ==="
rm -f namesorted/${SAMPLE}.ns.bam