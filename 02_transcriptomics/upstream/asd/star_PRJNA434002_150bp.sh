#!/bin/bash
#SBATCH --job-name=star_434002_150
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=logs/star_150bp_%j.out
#SBATCH --error=logs/star_150bp_%j.err

module purge
module load bear-apps/2024a
module load STAR/2.7.11b-GCC-13.3.0

GENOME=/rds/projects/e/elhamsak-dissertation-autism/RNA_seq/reference/GRCh38
cd /rds/projects/e/elhamsak-dissertation-autism/RNA_seq/PRJNA434002_bulk || exit 1
mkdir -p aligned logs

while read SAMPLE; do
  R1=trimmed/${SAMPLE}_1_val_1.fq.gz
  R2=trimmed/${SAMPLE}_2_val_2.fq.gz
  BAM=aligned/${SAMPLE}/Aligned.sortedByCoord.out.bam
  # Skip only if a REAL (non-empty) BAM exists
  if [ -f "$BAM" ] && [ $(stat -c%s "$BAM") -gt 1000000 ]; then
    echo "Skip (already aligned): $SAMPLE"
    continue
  fi
  if [ ! -f "$R1" ] || [ ! -f "$R2" ]; then
    echo "Missing trimmed: $SAMPLE"; continue
  fi
  mkdir -p aligned/${SAMPLE}
  echo "Aligning $SAMPLE ($(date))"
  STAR --runThreadN 16 --genomeDir $GENOME \
    --readFilesIn $R1 $R2 --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate --outSAMattributes NH HI AS NM \
    --outFileNamePrefix aligned/${SAMPLE}/ --quantMode GeneCounts --sjdbOverhang 99 \
    --outFilterMultimapNmax 20 --alignSJoverhangMin 8 \
    --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 \
    --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000
  echo "Done: $SAMPLE ($(date))"
done < metadata/keep_150bp.txt

echo "=== Complete: $(date) ==="
echo "Real BAMs now: $(for s in $(cat metadata/keep_150bp.txt); do b=aligned/$s/Aligned.sortedByCoord.out.bam; [ -f "$b" ] && [ $(stat -c%s "$b") -gt 1000000 ] && echo x; done | wc -l) / 41"