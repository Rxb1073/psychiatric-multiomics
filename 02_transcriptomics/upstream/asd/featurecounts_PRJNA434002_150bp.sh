#!/bin/bash
#SBATCH --job-name=fc_434002
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --output=logs/fc_150bp_%j.out
#SBATCH --error=logs/fc_150bp_%j.err

module purge
module load bear-apps/2022a
module load Subread/2.0.4-GCC-11.3.0

cd /rds/projects/e/elhamsak-dissertation-autism-2/RNA_seq/PRJNA434002_bulk || exit 1
GTF=/rds/projects/e/elhamsak-dissertation-autism/RNA_seq/reference/genome/Homo_sapiens.GRCh38.109.gtf
mkdir -p counts logs /rds/projects/e/elhamsak-dissertation-autism-2/tmp

# Build BAM list ONLY from the 41 clean 150bp samples
BAMS=""
while read s; do
  b="aligned/${s}/Aligned.sortedByCoord.out.bam"
  if [ -f "$b" ] && [ $(stat -c%s "$b") -gt 1000000 ]; then
    BAMS="$BAMS $b"
  fi
done < metadata/keep_150bp.txt

echo "Number of BAMs: $(echo $BAMS | wc -w) (expect 41)"

featureCounts -T 16 -p --countReadPairs -s 2 \
  --tmpDir /rds/projects/e/elhamsak-dissertation-autism-2/tmp \
  -a $GTF \
  -o counts/PRJNA434002_150bp_gene_counts.txt \
  $BAMS

echo "Done: $(date)"
