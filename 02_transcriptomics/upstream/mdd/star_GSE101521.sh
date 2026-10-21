#!/bin/bash
#SBATCH --job-name=star_GSE101521
#SBATCH --account=elhamsak-2-dissertat
#SBATCH --qos=bbdefault
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --array=1-59%6
#SBATCH --output=logs/star_%A_%a.out
#SBATCH --error=logs/star_%A_%a.err
module purge; module load bear-apps/2024a
module load STAR/2.7.11b-GCC-13.3.0
module load SAMtools/1.21-GCC-12.3.0
DS=/rds/projects/e/elhamsak-dissertation-autism-2/MDD/transcriptomics/GSE101521
GENOME=/rds/projects/e/elhamsak-dissertation-autism/RNA_seq/reference/GRCh38
cd $DS || exit 1
mkdir -p aligned logs
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" metadata/samples.txt)
[ -z "$SAMPLE" ] && { echo "no sample at index $SLURM_ARRAY_TASK_ID"; exit 1; }
FINAL=aligned/$SAMPLE/Aligned.sortedByCoord.out.bam
if [ -f "$FINAL" ] && [ $(stat -c%s "$FINAL") -gt 1000000 ]; then echo "skip $SAMPLE"; exit 0; fi
rm -rf aligned/$SAMPLE; mkdir -p aligned/$SAMPLE
TMP=/rds/projects/e/elhamsak-dissertation-autism-2/MDD/tmp/star_${SLURM_ARRAY_TASK_ID}_$$; rm -rf $TMP
# STAR -> UNSORTED bam (bypasses STAR's failing internal sorter)
STAR --runThreadN 16 --genomeDir $GENOME \
  --readFilesIn trimmed/${SAMPLE}_1_val_1.fq.gz trimmed/${SAMPLE}_2_val_2.fq.gz \
  --readFilesCommand zcat --outFileNamePrefix aligned/$SAMPLE/ \
  --outSAMtype BAM Unsorted --outSAMattributes NH HI AS NM \
  --quantMode GeneCounts --sjdbOverhang 99 \
  --outTmpDir $TMP \
  --outFilterMultimapNmax 20 --alignSJoverhangMin 8 \
  --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 \
  --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000
rm -rf $TMP
# samtools coordinate-sort -> same artifact as ASD (Aligned.sortedByCoord.out.bam)
samtools sort -@ 16 -m 3G -T $DS/aligned/$SAMPLE/cs_$SAMPLE \
  -o aligned/$SAMPLE/Aligned.sortedByCoord.out.bam aligned/$SAMPLE/Aligned.out.bam
# verify sorted bam is real, then remove the unsorted intermediate
if [ $(stat -c%s "aligned/$SAMPLE/Aligned.sortedByCoord.out.bam") -gt 1000000 ]; then
  rm -f aligned/$SAMPLE/Aligned.out.bam
else
  echo "SORT FAILED for $SAMPLE"; exit 1
fi
