#!/bin/bash
#SBATCH --job-name=GSE57369_ChIP
#SBATCH --account=elhamsak-dissertation-autism
#SBATCH --qos=bbdefault
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=/rds/projects/e/elhamsak-dissertation-autism/Epigenomics/logs/GSE57369_%j.out
#SBATCH --error=/rds/projects/e/elhamsak-dissertation-autism/Epigenomics/logs/GSE57369_%j.err
 
# =============================================================================
#  JOB 1: GSE57369 — CHD8 ChIP-seq
#  iPSC-derived Neural Progenitor Cells (hg19, Single-end)
#  Output: epigenomics-CHD8/GSE57369/
# =============================================================================
 
set -euo pipefail
THREADS=${SLURM_CPUS_PER_TASK}
 
echo "======================================"
echo " Job: ${SLURM_JOB_ID} | GSE57369"
echo " Start: $(date)"
echo "======================================"
 
# ── Paths ──────────────────────────────────────────────────────────────────
REFS="/rds/projects/e/elhamsak-dissertation-autism/Epigenomics/refs"
DIR="/rds/projects/e/elhamsak-dissertation-autism/Epigenomics/epigenomics-CHD8/GSE57369"
GENOME_INDEX="${REFS}/bowtie2/hg19/hg19"
BLACKLIST="${REFS}/blacklists/hg19-blacklist.v2.bed"
mkdir -p "${DIR}"/{fastq,fastqc,trimmed,bam,peaks,bigwig,logs,qc}
 
# ── Samples ────────────────────────────────────────────────────────────────
CHIP_SAMPLES=(SRR1273688 SRR1273689 SRR1273695 SRR1273691 SRR1273692)
INPUT_SAMPLE="SRR1273693"
ALL_SAMPLES=("${CHIP_SAMPLES[@]}" "${INPUT_SAMPLE}")
 
# ==========================================================================
#  BLOCK 1: Download → Filter (bear-apps/2023a)
# ==========================================================================
module purge
module load bear-apps/2023a/live
module load SRA-Toolkit/3.0.10-gompi-2023a
module load FastQC/0.11.9-Java-11
module load Trim_Galore/0.6.10-GCCcore-12.3.0
module load Bowtie2/2.5.4-GCC-12.3.0
module load SAMtools/1.21-GCC-12.3.0
module load picard/2.25.1-Java-11
module load BEDTools/2.31.1-GCC-12.3.0
 
echo "[GSE57369] STEP 1: Downloading..."
for SRR in "${ALL_SAMPLES[@]}"; do
    [[ -f "${DIR}/fastq/.${SRR}.done" ]] && { echo "  Skip: ${SRR}"; continue; }
    prefetch "${SRR}" --output-directory "${DIR}/fastq" \
        2>> "${DIR}/logs/${SRR}_prefetch.log"
    fasterq-dump "${DIR}/fastq/${SRR}/${SRR}.sra" \
        --outdir "${DIR}/fastq" --threads "${THREADS}" \
        --progress 2>> "${DIR}/logs/${SRR}_fasterq.log"
    gzip "${DIR}/fastq/${SRR}.fastq"
    touch "${DIR}/fastq/.${SRR}.done"
    echo "  Done: ${SRR}"
done
 
echo "[GSE57369] STEP 2: FastQC..."
fastqc "${DIR}/fastq/"*.fastq.gz \
    --outdir "${DIR}/fastqc" --threads "${THREADS}" \
    2>> "${DIR}/logs/fastqc.log"
 
echo "[GSE57369] STEP 3: Trimming..."
for SRR in "${ALL_SAMPLES[@]}"; do
    trim_galore --cores "${THREADS}" --fastqc \
        --output_dir "${DIR}/trimmed" \
        "${DIR}/fastq/${SRR}.fastq.gz" \
        2>> "${DIR}/logs/${SRR}_trim.log"
done
 
echo "[GSE57369] STEP 4: Aligning to hg19..."
for SRR in "${ALL_SAMPLES[@]}"; do
    bowtie2 -x "${GENOME_INDEX}" \
        -U "${DIR}/trimmed/${SRR}_trimmed.fq.gz" \
        --no-unal -p "${THREADS}" \
        2>> "${DIR}/logs/${SRR}_bowtie2.log" \
    | samtools sort -@ "${THREADS}" -o "${DIR}/bam/${SRR}.sorted.bam"
    samtools index "${DIR}/bam/${SRR}.sorted.bam"
    samtools flagstat "${DIR}/bam/${SRR}.sorted.bam" \
        > "${DIR}/qc/${SRR}_flagstat.txt"
done
 
echo "[GSE57369] STEP 5: Deduplication..."
for SRR in "${ALL_SAMPLES[@]}"; do
    java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
        INPUT="${DIR}/bam/${SRR}.sorted.bam" \
        OUTPUT="${DIR}/bam/${SRR}.dedup.bam" \
        METRICS_FILE="${DIR}/qc/${SRR}_dup_metrics.txt" \
        REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT \
        2>> "${DIR}/logs/${SRR}_picard.log"
    samtools index "${DIR}/bam/${SRR}.dedup.bam"
done
 
echo "[GSE57369] STEP 6: Filtering..."
for SRR in "${ALL_SAMPLES[@]}"; do
    samtools view -@ "${THREADS}" -b -q 30 \
        "${DIR}/bam/${SRR}.dedup.bam" \
    | bedtools intersect -v -abam stdin -b "${BLACKLIST}" \
    > "${DIR}/bam/${SRR}.filtered.bam"
    samtools index "${DIR}/bam/${SRR}.filtered.bam"
done
 
# ==========================================================================
#  BLOCK 2: Peak calling + BigWigs + QC (bear-apps/2022a)
# ==========================================================================
module purge
module load bear-apps/2022a/live
module load MACS2/2.2.9.1-foss-2022a
module load deepTools/3.5.2-foss-2022a
 
echo "[GSE57369] STEP 7: Peak calling..."
for SRR in "${CHIP_SAMPLES[@]}"; do
    macs2 callpeak \
        -t "${DIR}/bam/${SRR}.filtered.bam" \
        -c "${DIR}/bam/${INPUT_SAMPLE}.filtered.bam" \
        -f BAM -g hs \
        --name "GSE57369_${SRR}" \
        --outdir "${DIR}/peaks" \
        -q 0.05 --keep-dup auto \
        2>> "${DIR}/logs/${SRR}_macs2.log"
done
 
echo "[GSE57369] STEP 8: BigWigs..."
for SRR in "${CHIP_SAMPLES[@]}"; do
    bamCompare \
        -b1 "${DIR}/bam/${SRR}.filtered.bam" \
        -b2 "${DIR}/bam/${INPUT_SAMPLE}.filtered.bam" \
        --binSize 10 --normalizeUsing RPKM \
        --operation log2 --scaleFactorsMethod None \
        -p "${THREADS}" \
        -o "${DIR}/bigwig/${SRR}_log2FC.bw" \
        2>> "${DIR}/logs/${SRR}_bamcompare.log"
done
 
echo "[GSE57369] STEP 9: QC..."
BAMS=()
for SRR in "${ALL_SAMPLES[@]}"; do
    BAMS+=("${DIR}/bam/${SRR}.filtered.bam")
done
 
plotFingerprint -b "${BAMS[@]}" -p "${THREADS}" \
    --plotFile "${DIR}/qc/fingerprint.png" \
    --outRawCounts "${DIR}/qc/fingerprint_counts.tab" \
    2>> "${DIR}/logs/fingerprint.log"
 
multiBamSummary bins -b "${BAMS[@]}" -p "${THREADS}" \
    -o "${DIR}/qc/bam_summary.npz" \
    2>> "${DIR}/logs/multibamsummary.log"
 
plotCorrelation -in "${DIR}/qc/bam_summary.npz" \
    --corMethod pearson --skipZeros \
    --whatToPlot heatmap --colorMap RdYlBu --plotNumbers \
    -o "${DIR}/qc/correlation_heatmap.png" \
    2>> "${DIR}/logs/plotcorrelation.log"
 
echo "======================================"
echo " GSE57369 COMPLETE: $(date)"
echo " Peaks:   ${DIR}/peaks/"
echo " BigWigs: ${DIR}/bigwig/"
echo " QC:      ${DIR}/qc/"
echo "======================================"
 