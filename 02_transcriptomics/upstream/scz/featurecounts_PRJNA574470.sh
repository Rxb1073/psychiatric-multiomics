#!/bin/bash
#SBATCH --job-name=fc_PRJNA574470
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --output=logs/fc_%j.out
#SBATCH --error=logs/fc_%j.err
module purge
module load bear-apps/2022a
module load Subread/2.0.4-GCC-11.3.0
cd /rds/projects/e/elhamsak-dissertation-neuropsychiatric/RNA_seq/PRJNA574470 || exit 1
mkdir -p counts logs tmp

# PE branch: -p --countReadPairs, -s 2
echo "=== PE featureCounts: $(date) ==="
featureCounts -T 16 -p --countReadPairs -s 2 --primary   --tmpDir /rds/projects/e/elhamsak-dissertation-neuropsychiatric2/RNA_seq/PRJNA574470/tmp   -a /rds/projects/e/elhamsak-dissertation-autism/RNA_seq/reference/genome/Homo_sapiens.GRCh38.109.gtf   -o counts/PRJNA574470_PE_gene_counts.txt   $(cat metadata/pe_bams.txt)

# SE branch: NO -p, -s 2
echo "=== SE featureCounts: $(date) ==="
featureCounts -T 16 -s 2 --primary   --tmpDir /rds/projects/e/elhamsak-dissertation-neuropsychiatric2/RNA_seq/PRJNA574470/tmp   -a /rds/projects/e/elhamsak-dissertation-autism/RNA_seq/reference/genome/Homo_sapiens.GRCh38.109.gtf   -o counts/PRJNA574470_SE_gene_counts.txt   $(cat metadata/se_bams.txt)
echo "=== both branches done: $(date) ==="
