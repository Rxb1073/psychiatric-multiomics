#!/bin/bash
#SBATCH --job-name=star_PRJNA434002
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=logs/star_%j.out
#SBATCH --error=logs/star_%j.err

module purge
module load bear-apps/2024a
module load STAR/2.7.11b-GCC-13.3.0

GENOME=/rds/projects/e/elhamsak-dissertation-autism/RNA_seq/reference/GRCh38
cd /rds/projects/e/elhamsak-dissertation-autism/RNA_seq/PRJNA434002_bulk || exit 1
mkdir -p aligned logs

echo "=== STAR Alignment Started: $(date) ==="

while read SAMPLE; do
  R1=trimmed/${SAMPLE}_1_val_1.fq.gz
  R2=trimmed/${SAMPLE}_2_val_2.fq.gz
  if [ ! -f "$R1" ] || [ ! -f "$R2" ]; then
    echo "Missing trimmed files for $SAMPLE - skipping"
    continue
  fi
  if [ ! -f "aligned/${SAMPLE}/Aligned.sortedByCoord.out.bam" ]; then
    mkdir -p aligned/${SAMPLE}
    STAR \
      --runThreadN 16 \
      --genomeDir $GENOME \
      --readFilesIn $R1 $R2 \
      --readFilesCommand zcat \
      --outSAMtype BAM SortedByCoordinate \
      --outSAMattributes NH HI AS NM \
      --outFileNamePrefix aligned/${SAMPLE}/ \
      --quantMode GeneCounts \
      --sjdbOverhang 99 \
      --outFilterMultimapNmax 20 \
      --alignSJoverhangMin 8 \
      --outFilterMismatchNmax 999 \
      --outFilterMismatchNoverReadLmax 0.04 \
      --alignIntronMin 20 \
      --alignIntronMax 1000000 \
      --alignMatesGapMax 1000000
    echo "Done: $SAMPLE"
  else
    echo "Skipping $SAMPLE - already aligned"
  fi
done < metadata/keep_samples.txt

echo "=== STAR Complete: $(date) ==="
echo "Aligned: $(ls aligned/*/Aligned.sortedByCoord.out.bam 2>/dev/null | wc -l)"