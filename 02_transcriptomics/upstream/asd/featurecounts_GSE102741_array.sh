#!/bin/bash
#SBATCH --job-name=fc_GSE102741_array
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --array=0-3
#SBATCH --output=logs/fc_array_%A_%a.out
#SBATCH --error=logs/fc_array_%A_%a.err

module purge
module load bear-apps/2022a
module load Subread/2.0.4-GCC-11.3.0

GTF=/rds/projects/e/elhamsak-dissertation-autism/RNA_seq/reference/genome/Homo_sapiens.GRCh38.109.gtf
cd /rds/projects/e/elhamsak-dissertation-autism/RNA_seq/GSE102741 || exit 1
mkdir -p counts logs /rds/projects/e/elhamsak-dissertation-autism/tmp/fc_batch_${SLURM_ARRAY_TASK_ID}

# Get all BAM files into an array, split into 4 groups by array task ID
mapfile -t ALL_BAMS < <(ls aligned/*/Aligned.sortedByCoord.out.bam | sort)
TOTAL=${#ALL_BAMS[@]}

BATCH_BAMS=()
for ((i=SLURM_ARRAY_TASK_ID; i<TOTAL; i+=4)); do
  BATCH_BAMS+=("${ALL_BAMS[$i]}")
done

echo "=== Batch $SLURM_ARRAY_TASK_ID: ${#BATCH_BAMS[@]} files ==="
printf '%s\n' "${BATCH_BAMS[@]}"

featureCounts \
  -T 1 \
  -p \
  --countReadPairs \
  -s 2 \
  --tmpDir /rds/projects/e/elhamsak-dissertation-autism/tmp/fc_batch_${SLURM_ARRAY_TASK_ID} \
  -a $GTF \
  -o counts/GSE102741_gene_counts_batch${SLURM_ARRAY_TASK_ID}.txt \
  "${BATCH_BAMS[@]}"

echo "=== Batch $SLURM_ARRAY_TASK_ID Complete: $(date) ==="