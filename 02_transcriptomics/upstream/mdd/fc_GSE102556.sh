#!/bin/bash
#SBATCH --job-name=fc_GSE102556
#SBATCH --account=elhamsak-2-dissertat
#SBATCH --qos=bbdefault
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --array=1-48%8
#SBATCH --output=/rds/projects/e/elhamsak-dissertation-autism-2/MDD/transcriptomics/GSE102556/logs/fc_%A_%a.out
#SBATCH --error=/rds/projects/e/elhamsak-dissertation-autism-2/MDD/transcriptomics/GSE102556/logs/fc_%A_%a.err
module purge; module load bear-apps/2022a
module load SAMtools/1.16.1-GCC-11.3.0
module load Subread/2.0.4-GCC-11.3.0
DS=/rds/projects/e/elhamsak-dissertation-autism-2/MDD/transcriptomics/GSE102556
GTF=/rds/projects/e/elhamsak-dissertation-autism/RNA_seq/reference/genome/Homo_sapiens.GRCh38.109.gtf
STRAND=1
cd $DS || exit 1
mkdir -p counts/perbam namesorted logs
BAM=$(sed -n "${SLURM_ARRAY_TASK_ID}p" metadata/bam_list.txt)
[ -z "$BAM" ] && { echo "no BAM at $SLURM_ARRAY_TASK_ID"; exit 0; }
SAMPLE=$(basename $(dirname "$BAM"))
OUT=counts/perbam/${SAMPLE}.txt
if [ -f "$OUT" ] && [ $(stat -c%s "$OUT") -gt 1000 ]; then echo "skip $SAMPLE"; exit 0; fi
TMP=/rds/projects/e/elhamsak-dissertation-autism-2/MDD/tmp/fc_${SLURM_ARRAY_TASK_ID}_$$; mkdir -p $TMP
samtools quickcheck "$BAM" || { echo "BAD BAM $SAMPLE"; exit 1; }
samtools sort -n -@ 8 -m 4G -T $TMP/ns -o namesorted/${SAMPLE}.ns.bam "$BAM"
featureCounts -T 8 -p --countReadPairs -s $STRAND --tmpDir $TMP \
  -a $GTF -o "$OUT" namesorted/${SAMPLE}.ns.bam
rm -f namesorted/${SAMPLE}.ns.bam
rm -rf $TMP
G=$(grep -vc '^#\|^Geneid' "$OUT")
echo "$SAMPLE: $G genes"
