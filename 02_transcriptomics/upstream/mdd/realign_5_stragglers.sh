#!/bin/bash
#SBATCH --job-name=realign5
#SBATCH --account=elhamsak-2-dissertat
#SBATCH --qos=bbdefault
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --array=1-5
#SBATCH --output=/rds/projects/e/elhamsak-dissertation-autism-2/MDD/transcriptomics/logs/realign5_%A_%a.out
#SBATCH --error=/rds/projects/e/elhamsak-dissertation-autism-2/MDD/transcriptomics/logs/realign5_%A_%a.err
set -uo pipefail
A2=/rds/projects/e/elhamsak-dissertation-autism-2/MDD/transcriptomics
GENOME=/rds/projects/e/elhamsak-dissertation-autism/RNA_seq/reference/GRCh38

# map array task -> "dataset sample"
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $A2/logs/straggler_map.txt)
D=$(echo $LINE | awk '{print $1}')
SAMPLE=$(echo $LINE | awk '{print $2}')
DS=$A2/$D
echo "=== task $SLURM_ARRAY_TASK_ID : $D $SAMPLE ==="

cd $DS || exit 1
rm -rf aligned/$SAMPLE; mkdir -p aligned/$SAMPLE

# --- STAR (2024a) -> unsorted BAM ---
module purge; module load bear-apps/2024a; module load STAR/2.7.11b-GCC-13.3.0
TMP=$A2/../tmp/ra_${SAMPLE}_$$; rm -rf $TMP
STAR --runThreadN 16 --genomeDir $GENOME \
  --readFilesIn trimmed/${SAMPLE}_1_val_1.fq.gz trimmed/${SAMPLE}_2_val_2.fq.gz \
  --readFilesCommand zcat --outFileNamePrefix aligned/$SAMPLE/ \
  --outSAMtype BAM Unsorted --outSAMattributes NH HI AS NM \
  --quantMode GeneCounts --sjdbOverhang 99 --outTmpDir $TMP \
  --outFilterMultimapNmax 20 --alignSJoverhangMin 8 \
  --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 \
  --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000
rm -rf $TMP
if [ ! -s aligned/$SAMPLE/Aligned.out.bam ]; then echo "STAR FAILED $SAMPLE"; exit 1; fi

# --- samtools (2022a) -> coordinate sort ---
module purge; module load bear-apps/2022a; module load SAMtools/1.16.1-GCC-11.3.0
TMP2=$A2/../tmp/ras_${SAMPLE}_$$; mkdir -p $TMP2
samtools sort -@ 16 -m 3G -T $TMP2/s \
  -o aligned/$SAMPLE/Aligned.sortedByCoord.out.bam aligned/$SAMPLE/Aligned.out.bam
# FULL read verification (not just quickcheck)
if samtools view -c aligned/$SAMPLE/Aligned.sortedByCoord.out.bam >/dev/null 2>&1; then
  rm -f aligned/$SAMPLE/Aligned.out.bam
  echo "SORT+VERIFY OK $SAMPLE"
else
  echo "SORT/VERIFY FAILED $SAMPLE"; exit 1
fi
rm -rf $TMP2
