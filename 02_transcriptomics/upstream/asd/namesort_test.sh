#!/bin/bash
#SBATCH --job-name=namesort_test
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --output=logs/namesort_%j.out
#SBATCH --error=logs/namesort_%j.err
module purge
module load bear-apps/2022a
module load SAMtools/1.16.1-GCC-11.3.0
module load Subread/2.0.4-GCC-11.3.0
GSE=/rds/projects/e/elhamsak-dissertation-autism/RNA_seq/GSE102741
GTF=/rds/projects/e/elhamsak-dissertation-autism/RNA_seq/reference/genome/Homo_sapiens.GRCh38.109.gtf
cd $GSE || exit 1
mkdir -p counts/perbam logs namesorted
TMP=/rds/projects/e/elhamsak-dissertation-autism/tmp/ns_test
mkdir -p $TMP

SAMPLE=SRR5938441   # the 36GB stuck BAM (task 22)
BAM=$GSE/aligned/$SAMPLE/Aligned.sortedByCoord.out.bam

echo "=== Name-sorting $SAMPLE ($(date)) BAM size: $(ls -lh $BAM | awk '{print $5}') ==="
samtools sort -n -@ 8 -m 4G -T $TMP/ns_$SAMPLE \
  -o namesorted/${SAMPLE}.namesorted.bam "$BAM"
echo "=== Name-sort done ($(date)). Size: $(ls -lh namesorted/${SAMPLE}.namesorted.bam | awk '{print $5}') ==="

echo "=== featureCounts on name-sorted BAM ($(date)) ==="
featureCounts -T 8 -p --countReadPairs -s 2 \
  --tmpDir $TMP \
  -a $GTF \
  -o counts/perbam/${SAMPLE}.txt \
  namesorted/${SAMPLE}.namesorted.bam
echo "=== Done ($(date)) ==="
rm -f namesorted/${SAMPLE}.namesorted.bam   # cleanup the big intermediate