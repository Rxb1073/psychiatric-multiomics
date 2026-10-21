#!/bin/bash
#SBATCH --job-name=fc_GSE102741
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --output=logs/fc_%j.out
#SBATCH --error=logs/fc_%j.err

module purge
module load bear-apps/2022a
module load Subread/2.0.4-GCC-11.3.0

GTF=/rds/projects/e/elhamsak-dissertation-autism/RNA_seq/reference/genome/Homo_sapiens.GRCh38.109.gtf
cd /rds/projects/e/elhamsak-dissertation-autism/RNA_seq/GSE102741 || exit 1
mkdir -p counts logs

featureCounts -T 16 -p --countReadPairs -s 2 \
  -a $GTF \
  -o counts/GSE102741_gene_counts.txt \
  aligned/*/Aligned.sortedByCoord.out.bam

echo "Done: $(date)"