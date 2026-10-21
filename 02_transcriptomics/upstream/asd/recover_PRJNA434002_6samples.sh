#!/bin/bash
#SBATCH --job-name=recover_prjna
#SBATCH --array=0-5
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --output=logs/recover_%A_%a.out
#SBATCH --error=logs/recover_%A_%a.err
# ---------------------------------------------------------------------------
# Recover the 6 PRJNA434002 samples whose paired trimming produced empty
# output due to SRA read-name formatting ("/1 /2" suffix breaking pairing).
# Reads are intact and in-sync, so: repair -> Trim Galore -> STAR.
# Nothing is re-downloaded. Original raw files are read-only inputs.
# ---------------------------------------------------------------------------
set -euo pipefail

DS=/rds/projects/e/elhamsak-dissertation-autism-2/RNA_seq/PRJNA434002_bulk
GENOME=/rds/projects/e/elhamsak-dissertation-autism/RNA_seq/reference/GRCh38
GTF=/rds/projects/e/elhamsak-dissertation-autism/RNA_seq/reference/genome/Homo_sapiens.GRCh38.109.gtf
cd "$DS"

SAMPLES=(SRR9262917 SRR9262919 SRR9262920 SRR9262921 SRR9262923 SRR9262927)
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

TMP=/rds/projects/e/elhamsak-dissertation-autism-2/tmp/recover_${SAMPLE}_$$
mkdir -p "$TMP" repaired trimmed aligned/"$SAMPLE" logs
trap 'rm -rf "$TMP"' EXIT

echo "=== [$SAMPLE] START $(date) ==="

# 1) REPAIR: re-pair mates on base read-name, normalising the /1 /2 suffix ----
module purge
module load bear-apps/2024a
module load BBMap/39.19-GCC-13.3.0
echo "--- repair.sh ---"
repair.sh -Xmx40g \
  in1=fastq_files/${SAMPLE}_1.fastq.gz \
  in2=fastq_files/${SAMPLE}_2.fastq.gz \
  out1=repaired/${SAMPLE}_1.fastq.gz \
  out2=repaired/${SAMPLE}_2.fastq.gz \
  repair tossbrokenreads

# 2) TRIM GALORE: identical parameters to the rest of PRJNA434002 -------------
module purge
module load bear-apps/2023a
module load Trim_Galore/0.6.10-GCCcore-12.3.0
echo "--- trim_galore ---"
trim_galore \
  --quality 20 \
  --length 36 \
  --paired \
  --cores 4 \
  --output_dir trimmed \
  repaired/${SAMPLE}_1.fastq.gz repaired/${SAMPLE}_2.fastq.gz

# sanity check: the val files must now be non-empty
R1V=trimmed/${SAMPLE}_1_val_1.fq.gz
R2V=trimmed/${SAMPLE}_2_val_2.fq.gz
if [ "$(stat -c%s "$R1V")" -lt 1000 ] || [ "$(stat -c%s "$R2V")" -lt 1000 ]; then
  echo "!!! [$SAMPLE] trimmed output still empty — stopping before STAR" >&2
  exit 1
fi

# 3) STAR: same params + reference as the rest of the cohort (-s 2 confirmed) --
module purge
module load bear-apps/2024a
module load STAR/2.7.11b-GCC-13.3.0
echo "--- STAR ---"
STAR --runThreadN 16 --genomeDir "$GENOME" \
  --readFilesIn "$R1V" "$R2V" --readFilesCommand zcat \
  --outFileNamePrefix aligned/${SAMPLE}/ \
  --outSAMtype BAM SortedByCoordinate --outSAMattributes NH HI AS NM \
  --quantMode GeneCounts --sjdbOverhang 99 \
  --outTmpDir "$TMP"/star \
  --outFilterMultimapNmax 20 --alignSJoverhangMin 8 \
  --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 \
  --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000

echo "=== [$SAMPLE] DONE $(date) ==="
grep 'Uniquely mapped reads %' aligned/${SAMPLE}/Log.final.out