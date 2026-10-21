#!/bin/bash
#SBATCH --job-name=star_PRJNA314463
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=40G
#SBATCH --array=1-82%20
#SBATCH --output=/rds/projects/e/elhamsak-dissertation-neuropsychiatric2/RNA_seq/PRJNA314463/logs/star_%A_%a.out
#SBATCH --error=/rds/projects/e/elhamsak-dissertation-neuropsychiatric2/RNA_seq/PRJNA314463/logs/star_%A_%a.err
module purge
module load bear-apps/2024a
module load STAR/2.7.11b-GCC-13.3.0
SRC=/rds/projects/e/elhamsak-dissertation-neuropsychiatric/RNA_seq/PRJNA314463
OUT_VOL=/rds/projects/e/elhamsak-dissertation-neuropsychiatric2/RNA_seq/PRJNA314463
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $SRC/metadata/star_samples.txt)
OUT=$OUT_VOL/aligned/${SAMPLE}
if [ -f "$OUT/Aligned.sortedByCoord.out.bam" ] && [ $(stat -c%s "$OUT/Aligned.sortedByCoord.out.bam") -gt 1000000 ]; then
  echo "Skip $SAMPLE"; exit 0; fi
rm -rf $OUT && mkdir -p $OUT
TMP=$OUT_VOL/tmp/${SAMPLE}_${SLURM_ARRAY_TASK_ID}; rm -rf $TMP
echo "=== STAR $SAMPLE start: $(date) ==="
STAR --runThreadN 8 --genomeDir /rds/projects/e/elhamsak-dissertation-autism/RNA_seq/reference/GRCh38   --readFilesIn $SRC/trimmed/${SAMPLE}_1_val_1.fq.gz $SRC/trimmed/${SAMPLE}_2_val_2.fq.gz   --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate   --outSAMattributes NH HI AS NM --outFileNamePrefix $OUT/   --outTmpDir $TMP   --quantMode GeneCounts --sjdbOverhang 99   --outFilterMultimapNmax 20 --alignSJoverhangMin 8   --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04   --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000
echo "=== $SAMPLE done: $(date) ==="
