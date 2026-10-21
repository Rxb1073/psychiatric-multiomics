#!/bin/bash
#SBATCH --job-name=fc_stragglers
#SBATCH --account=elhamsak-2-dissertat
#SBATCH --qos=bbdefault
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --output=/rds/projects/e/elhamsak-dissertation-autism-2/MDD/transcriptomics/logs/fc_stragglers_%j.out
#SBATCH --error=/rds/projects/e/elhamsak-dissertation-autism-2/MDD/transcriptomics/logs/fc_stragglers_%j.err
module purge; module load bear-apps/2022a
module load SAMtools/1.16.1-GCC-11.3.0
module load Subread/2.0.4-GCC-11.3.0
GTF=/rds/projects/e/elhamsak-dissertation-autism/RNA_seq/reference/genome/Homo_sapiens.GRCh38.109.gtf
A2=/rds/projects/e/elhamsak-dissertation-autism-2/MDD/transcriptomics
run_one () {
  local D=$1 STRAND=$2 SAMPLE=$3
  local DS=$A2/$D
  local BAM=$DS/aligned/$SAMPLE/Aligned.sortedByCoord.out.bam
  local OUT=$DS/counts/perbam/${SAMPLE}.txt
  local TMP=$A2/../tmp/fcs_${SAMPLE}_$$; mkdir -p $TMP
  echo ">>> $D $SAMPLE strand=$STRAND"
  samtools view -c "$BAM" >/dev/null 2>&1 || { echo "BAD BAM $SAMPLE"; return 1; }
  samtools sort -n -@ 8 -m 4G -T $TMP/ns -o $TMP/${SAMPLE}.ns.bam "$BAM"
  featureCounts -T 8 -p --countReadPairs -s $STRAND --tmpDir $TMP \
    -a $GTF -o "$OUT" $TMP/${SAMPLE}.ns.bam
  echo "    $SAMPLE -> $(grep -vc '^#\|^Geneid' "$OUT") genes"
  rm -rf $TMP
}
run_one GSE102556 1 SRR5961874
run_one GSE102556 1 SRR5961878
run_one GSE101521 2 SRR5831972
run_one GSE101521 2 SRR5831974
run_one GSE101521 2 SRR5831975
echo "=== stragglers done ==="
