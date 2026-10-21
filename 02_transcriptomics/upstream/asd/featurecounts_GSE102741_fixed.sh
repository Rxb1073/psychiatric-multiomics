#!/bin/bash
#SBATCH --job-name=fc_GSE102741
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --output=logs/fc_fixed_%j.out
#SBATCH --error=logs/fc_fixed_%j.err

module purge
module load bear-apps/2022a
module load Subread/2.0.4-GCC-11.3.0

cd /rds/projects/e/elhamsak-dissertation-autism/RNA_seq/GSE102741 || exit 1
GTF=/rds/projects/e/elhamsak-dissertation-autism/RNA_seq/reference/genome/Homo_sapiens.GRCh38.109.gtf
mkdir -p counts logs /rds/projects/e/elhamsak-dissertation-autism/tmp

# Build BAM list from all real BAMs
BAMS=$(find aligned -name "Aligned.sortedByCoord.out.bam" -size +1M | sort)
echo "Number of BAMs: $(echo "$BAMS" | wc -l) (expect 52)"

featureCounts -T 16 -p --countReadPairs -s 2 \
  --tmpDir /rds/projects/e/elhamsak-dissertation-autism/tmp \
  -a $GTF \
  -o counts/GSE102741_gene_counts.txt \
  $BAMS

echo "Done: $(date)"
