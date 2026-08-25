#!/bin/bash
#SBATCH --job-name=GSE61492_ChIP
#SBATCH --account=elhamsak-dissertation-autism
#SBATCH --qos=bbdefault
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=/rds/projects/e/elhamsak-dissertation-autism/Epigenomics/logs/GSE61492_%j.out
#SBATCH --error=/rds/projects/e/elhamsak-dissertation-autism/Epigenomics/logs/GSE61492_%j.err

set -euo pipefail
THREADS=${SLURM_CPUS_PER_TASK}

echo "======================================"
echo " Job: ${SLURM_JOB_ID} | GSE61492"
echo " Start: $(date)"
echo "======================================"

module purge
module load bear-apps/2023a/live
module load SRA-Toolkit/3.0.10-gompi-2023a
module load FastQC/0.11.9-Java-11
module load Trim_Galore/0.6.10-GCCcore-12.3.0
module load Bowtie2/2.5.4-GCC-12.3.0
module load SAMtools/1.21-GCC-12.3.0
module load picard/2.25.1-Java-11
module load BEDTools/2.31.1-GCC-12.3.0

REFS="/rds/projects/e/elhamsak-dissertation-autism/Epigenomics/refs"
DIR="/rds/projects/e/elhamsak-dissertation-autism/Epigenomics/epigenomics-CHD8/GSE61492"
GENOME_INDEX="${REFS}/bowtie2/hg38/hg38"
BLACKLIST="${REFS}/blacklists/hg38-blacklist.v2.bed"
mkdir -p "${DIR}"/{fastq,fastqc,trimmed,bam,peaks,bigwig,logs,qc}

CHIP_SAMPLES=(SRR1576798 SRR1576799 SRR1576800 SRR1576801 SRR1576802 SRR1576803 SRR1576804)
INPUT_SAMPLE="SRR1576805"
ALL_SAMPLES=("${CHIP_SAMPLES[@]}" "${INPUT_SAMPLE}")

echo "[GSE61492] STEP 1: Downloading..."
for SRR in "${ALL_SAMPLES[@]}"; do
    [[ -f "${DIR}/fastq/.${SRR}.done" ]] && { echo "  Skip: ${SRR}"; continue; }
    prefetch "${SRR}" --output-directory "${DIR}/fastq" \
        2>> "${DIR}/logs/${SRR}_prefetch.log"
    # Use full path to sralite or sra file
    SRA_FILE=$(find "${DIR}/fastq/${SRR}/" \( -name "*.sra" -o -name "*.sralite" \) 2>/dev/null | head -1)
    fasterq-dump "${SRA_FILE}" \
        --outdir "${DIR}/fastq" --threads "${THREADS}" \
        --progress 2>> "${DIR}/logs/${SRR}_fasterq.log"
    # Handle both SE and PE output from fasterq-dump
    if [[ -f "${DIR}/fastq/${SRR}_1.fastq" ]]; then
        mv "${DIR}/fastq/${SRR}_1.fastq" "${DIR}/fastq/${SRR}.fastq"
        rm -f "${DIR}/fastq/${SRR}_2.fastq"
    fi
    gzip "${DIR}/fastq/${SRR}.fastq"
    touch "${DIR}/fastq/.${SRR}.done"
    echo "  Done: ${SRR}"
done

echo "[GSE61492] STEP 2: FastQC..."
fastqc "${DIR}/fastq/"*.fastq.gz \
    --outdir "${DIR}/fastqc" --threads "${THREADS}" \
    2>> "${DIR}/logs/fastqc.log"

echo "[GSE61492] STEP 3: Trimming..."
for SRR in "${ALL_SAMPLES[@]}"; do
    trim_galore --cores "${THREADS}" --fastqc \
        --output_dir "${DIR}/trimmed" \
        "${DIR}/fastq/${SRR}.fastq.gz" \
        2>> "${DIR}/logs/${SRR}_trim.log"
done

echo "[GSE61492] STEP 4: Aligning to hg38..."
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

echo "[GSE61492] STEP 5: Deduplication..."
for SRR in "${ALL_SAMPLES[@]}"; do
    java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
        INPUT="${DIR}/bam/${SRR}.sorted.bam" \
        OUTPUT="${DIR}/bam/${SRR}.dedup.bam" \
        METRICS_FILE="${DIR}/qc/${SRR}_dup_metrics.txt" \
        REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT \
        2>> "${DIR}/logs/${SRR}_picard.log"
    samtools index "${DIR}/bam/${SRR}.dedup.bam"
done

echo "[GSE61492] STEP 6: Filtering..."
for SRR in "${ALL_SAMPLES[@]}"; do
    samtools view -@ "${THREADS}" -b -q 30 \
        "${DIR}/bam/${SRR}.dedup.bam" \
    | bedtools intersect -v -abam stdin -b "${BLACKLIST}" \
    > "${DIR}/bam/${SRR}.filtered.bam"
    samtools index "${DIR}/bam/${SRR}.filtered.bam"
done

module purge
module load bear-apps/2022a/live
module load MACS2/2.2.9.1-foss-2022a
module load deepTools/3.5.2-foss-2022a

echo "[GSE61492] STEP 7: Peak calling..."
for SRR in "${CHIP_SAMPLES[@]}"; do
    macs2 callpeak \
        -t "${DIR}/bam/${SRR}.filtered.bam" \
        -c "${DIR}/bam/${INPUT_SAMPLE}.filtered.bam" \
        -f BAM -g hs --name "GSE61492_${SRR}" \
        --outdir "${DIR}/peaks" -q 0.05 --keep-dup auto \
        2>> "${DIR}/logs/${SRR}_macs2.log"
done

echo "[GSE61492] STEP 8: BigWigs..."
for SRR in "${CHIP_SAMPLES[@]}"; do
    bamCompare \
        -b1 "${DIR}/bam/${SRR}.filtered.bam" \
        -b2 "${DIR}/bam/${INPUT_SAMPLE}.filtered.bam" \
        --binSize 10 --normalizeUsing RPKM \
        --operation log2 --scaleFactorsMethod None \
        -p "${THREADS}" -o "${DIR}/bigwig/${SRR}_log2FC.bw" \
        2>> "${DIR}/logs/${SRR}_bamcompare.log"
done

echo "[GSE61492] STEP 9: QC..."
BAMS=()
for SRR in "${ALL_SAMPLES[@]}"; do BAMS+=("${DIR}/bam/${SRR}.filtered.bam"); done

plotFingerprint -b "${BAMS[@]}" -p "${THREADS}" \
    --plotFile "${DIR}/qc/fingerprint.png" \
    --outRawCounts "${DIR}/qc/fingerprint_counts.tab" \
    2>> "${DIR}/logs/fingerprint.log"

multiBamSummary bins -b "${BAMS[@]}" -p "${THREADS}" \
    -o "${DIR}/qc/bam_summary.npz" 2>> "${DIR}/logs/multibamsummary.log"

plotCorrelation -in "${DIR}/qc/bam_summary.npz" \
    --corMethod pearson --skipZeros \
    --whatToPlot heatmap --colorMap RdYlBu --plotNumbers \
    -o "${DIR}/qc/correlation_heatmap.png" \
    2>> "${DIR}/logs/plotcorrelation.log"

echo "======================================"
echo " GSE61492 COMPLETE: $(date)"
echo "======================================"
