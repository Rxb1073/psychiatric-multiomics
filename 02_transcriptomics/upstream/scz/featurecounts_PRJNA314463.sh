#!/bin/bash
#SBATCH --job-name=fc_PRJNA314463
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --output=logs/fc_%j.out
#SBATCH --error=logs/fc_%j.err
module purge
module load bear-apps/2022a
module load Subread/2.0.4-GCC-11.3.0
cd /rds/projects/e/elhamsak-dissertation-neuropsychiatric/RNA_seq/PRJNA314463 || exit 1
mkdir -p counts logs tmp
featureCounts -T 16 -p --countReadPairs -s 0 --primary   --tmpDir /rds/projects/e/elhamsak-dissertation-neuropsychiatric2/RNA_seq/PRJNA314463/tmp   -a /rds/projects/e/elhamsak-dissertation-autism/RNA_seq/reference/genome/Homo_sapiens.GRCh38.109.gtf   -o counts/PRJNA314463_gene_counts.txt   /rds/projects/e/elhamsak-dissertation-neuropsychiatric2/RNA_seq/PRJNA314463/aligned/*/Aligned.sortedByCoord.out.bam
echo "=== featureCounts done: $(date) ==="
