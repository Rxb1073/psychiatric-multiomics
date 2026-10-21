#!/bin/bash
#SBATCH --job-name=realign_fc_434
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=logs/realign_fc_%j.out
#SBATCH --error=logs/realign_fc_%j.err

cd /rds/projects/e/elhamsak-dissertation-autism-2/RNA_seq/PRJNA434002_bulk || exit 1
mkdir -p aligned logs counts
export TMPDIR=/rds/projects/e/elhamsak-dissertation-autism-2/tmp

GENOME=/rds/projects/e/elhamsak-dissertation-autism/RNA_seq/reference/GRCh38
GTF=/rds/projects/e/elhamsak-dissertation-autism/RNA_seq/reference/genome/Homo_sapiens.GRCh38.109.gtf

# --- Step 1: re-align SRR9292593 ---
module purge
module load bear-apps/2024a
module load STAR/2.7.11b-GCC-13.3.0

SAMPLE=SRR9292593
echo "Re-aligning $SAMPLE ($(date))"
mkdir -p aligned/${SAMPLE}
STAR --runThreadN 16 --genomeDir $GENOME \
  --readFilesIn trimmed/${SAMPLE}_1_val_1.fq.gz trimmed/${SAMPLE}_2_val_2.fq.gz \
  --readFilesCommand zcat \
  --outSAMtype BAM SortedByCoordinate --outSAMattributes NH HI AS NM \
  --outFileNamePrefix aligned/${SAMPLE}/ --quantMode GeneCounts --sjdbOverhang 99 \
  --outFilterMultimapNmax 20 --alignSJoverhangMin 8 \
  --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 \
  --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000
echo "Re-align done ($(date))"

# --- Step 2: featureCounts on all 41 ---
module purge
module load bear-apps/2022a
module load Subread/2.0.4-GCC-11.3.0

BAMS=""
n=0
while read s; do
  b="aligned/${s}/Aligned.sortedByCoord.out.bam"
  if [ -f "$b" ] && [ $(stat -c%s "$b") -gt 1000000 ]; then
    BAMS="$BAMS $b"; n=$((n+1))
  fi
done < metadata/keep_150bp.txt
echo "featureCounts on $n BAMs (expect 41)"

featureCounts -T 16 -p --countReadPairs -s 2 \
  --tmpDir $TMPDIR \
  -a $GTF \
  -o counts/PRJNA434002_150bp_gene_counts.txt \
  $BAMS

echo "=== Complete: $(date) ==="
