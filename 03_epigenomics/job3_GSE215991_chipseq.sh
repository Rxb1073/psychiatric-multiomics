#!/bin/bash
#SBATCH --job-name=GSE215991_ChIP
#SBATCH --account=elhamsak-4-dissertat
#SBATCH --qos=bbdefault
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output=/rds/projects/e/elhamsak-dissertation-neuropsychiatric2/Epigenomics/logs/GSE215991_%j.out
#SBATCH --error=/rds/projects/e/elhamsak-dissertation-neuropsychiatric2/Epigenomics/logs/GSE215991_%j.err

set -euo pipefail
THREADS=${SLURM_CPUS_PER_TASK}

echo "======================================"
echo " Job: ${SLURM_JOB_ID} | GSE215991"
echo " Start: $(date)"
echo "======================================"

REFS="/rds/projects/e/elhamsak-dissertation-autism/Epigenomics/refs"
DIR="/rds/projects/e/elhamsak-dissertation-autism/Epigenomics/epigenomics-ChIP-seq/GSE215991"
GENOME_INDEX="${REFS}/bowtie2/hg38/hg38"
BLACKLIST="${REFS}/blacklists/hg38-blacklist.v2.bed"
mkdir -p "${DIR}"/{fastq,fastqc,trimmed,bam,peaks,bigwig,logs,qc}

SCZ_SAMPLES=(SRR21948648 SRR21948649 SRR21948650 SRR21948651 SRR21948652 SRR21948653 SRR21948654 SRR21948655 SRR21948656 SRR21948657 SRR21948658 SRR21948659 SRR21948660 SRR21948661 SRR21948662)
NPC_SAMPLES=(SRR21948663 SRR21948664 SRR21948665 SRR21948666 SRR21948667 SRR21948668 SRR21948669 SRR21948670 SRR21948671 SRR21948672 SRR21948673 SRR21948674 SRR21948675 SRR21948676 SRR21948677)
ALL_SAMPLES=("${SCZ_SAMPLES[@]}" "${NPC_SAMPLES[@]}")

module purge
module load bear-apps/2023a/live
module load SRA-Toolkit/3.0.10-gompi-2023a
module load FastQC/0.11.9-Java-11
module load Trim_Galore/0.6.10-GCCcore-12.3.0
module load Bowtie2/2.5.4-GCC-12.3.0
module load SAMtools/1.21-GCC-12.3.0
module load picard/2.25.1-Java-11
module load BEDTools/2.31.1-GCC-12.3.0

echo "[GSE215991] STEP 1: Downloading 30 samples..."
for SRR in "${ALL_SAMPLES[@]}"; do
    [[ -f "${DIR}/fastq/.${SRR}.done" ]] && { echo "  Skip: ${SRR}"; continue; }
    prefetch "${SRR}" --output-directory "${DIR}/fastq" 2>> "${DIR}/logs/${SRR}_prefetch.log"
    SRA_FILE=$(find "${DIR}/fastq/${SRR}/" \( -name "*.sra" -o -name "*.sralite" \) 2>/dev/null | head -1)
    fasterq-dump "${SRA_FILE}" --outdir "${DIR}/fastq" --threads "${THREADS}" --progress 2>> "${DIR}/logs/${SRR}_fasterq.log"
    if [[ -f "${DIR}/fastq/${SRR}_1.fastq" ]]; then
        mv "${DIR}/fastq/${SRR}_1.fastq" "${DIR}/fastq/${SRR}.fastq"
        rm -f "${DIR}/fastq/${SRR}_2.fastq"
    fi
    gzip "${DIR}/fastq/${SRR}.fastq"
    touch "${DIR}/fastq/.${SRR}.done"
    echo "  Done: ${SRR}"
done

echo "[GSE215991] STEP 2: FastQC..."
fastqc "${DIR}/fastq/"*.fastq.gz --outdir "${DIR}/fastqc" --threads "${THREADS}" 2>> "${DIR}/logs/fastqc.log"

echo "[GSE215991] STEP 3: Trimming..."
for SRR in "${ALL_SAMPLES[@]}"; do
    trim_galore --cores "${THREADS}" --fastqc --output_dir "${DIR}/trimmed" "${DIR}/fastq/${SRR}.fastq.gz" 2>> "${DIR}/logs/${SRR}_trim.log"
done

echo "[GSE215991] STEP 4: Aligning..."
for SRR in "${ALL_SAMPLES[@]}"; do
    bowtie2 -x "${GENOME_INDEX}" -U "${DIR}/trimmed/${SRR}_trimmed.fq.gz" --no-unal -p "${THREADS}" 2>> "${DIR}/logs/${SRR}_bowtie2.log" | samtools sort -@ "${THREADS}" -o "${DIR}/bam/${SRR}.sorted.bam"
    samtools index "${DIR}/bam/${SRR}.sorted.bam"
    samtools flagstat "${DIR}/bam/${SRR}.sorted.bam" > "${DIR}/qc/${SRR}_flagstat.txt"
    echo "  Aligned: ${SRR}"
done

echo "[GSE215991] STEP 5: Deduplication..."
for SRR in "${ALL_SAMPLES[@]}"; do
    java -jar $EBROOTPICARD/picard.jar MarkDuplicates INPUT="${DIR}/bam/${SRR}.sorted.bam" OUTPUT="${DIR}/bam/${SRR}.dedup.bam" METRICS_FILE="${DIR}/qc/${SRR}_dup_metrics.txt" REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT 2>> "${DIR}/logs/${SRR}_picard.log"
    samtools index "${DIR}/bam/${SRR}.dedup.bam"
done

echo "[GSE215991] STEP 6: Filtering..."
for SRR in "${ALL_SAMPLES[@]}"; do
    samtools view -@ "${THREADS}" -b -q 30 "${DIR}/bam/${SRR}.dedup.bam" | bedtools intersect -v -abam stdin -b "${BLACKLIST}" > "${DIR}/bam/${SRR}.filtered.bam"
    samtools index "${DIR}/bam/${SRR}.filtered.bam"
done

module purge
module load bear-apps/2022a/live
module load MACS2/2.2.9.1-foss-2022a
module load deepTools/3.5.2-foss-2022a

echo "[GSE215991] STEP 7: Peak calling..."
SCZ_BAMS=(); for SRR in "${SCZ_SAMPLES[@]}"; do SCZ_BAMS+=("${DIR}/bam/${SRR}.filtered.bam"); done
NPC_BAMS=(); for SRR in "${NPC_SAMPLES[@]}"; do NPC_BAMS+=("${DIR}/bam/${SRR}.filtered.bam"); done

macs2 callpeak -t "${SCZ_BAMS[@]}" -c "${NPC_BAMS[@]}" -f BAM -g hs --name "GSE215991_SCZ_vs_NPC" --outdir "${DIR}/peaks" --broad --broad-cutoff 0.1 --keep-dup auto 2>> "${DIR}/logs/macs2_pooled.log"

echo "[GSE215991] STEP 8: BigWigs..."
module purge
module load bear-apps/2023a/live
module load SAMtools/1.21-GCC-12.3.0
samtools merge -@ "${THREADS}" -f "${DIR}/bam/SCZ_merged.bam" "${SCZ_BAMS[@]}"
samtools index "${DIR}/bam/SCZ_merged.bam"
samtools merge -@ "${THREADS}" -f "${DIR}/bam/NPC_merged.bam" "${NPC_BAMS[@]}"
samtools index "${DIR}/bam/NPC_merged.bam"

module purge
module load bear-apps/2022a/live
module load deepTools/3.5.2-foss-2022a
bamCompare -b1 "${DIR}/bam/SCZ_merged.bam" -b2 "${DIR}/bam/NPC_merged.bam" --binSize 10 --normalizeUsing RPKM --operation log2 --scaleFactorsMethod None -p "${THREADS}" -o "${DIR}/bigwig/SCZ_vs_NPC_log2FC.bw" 2>> "${DIR}/logs/bamcompare.log"

echo "[GSE215991] STEP 9: QC..."
plotFingerprint -b "${DIR}/bam/SCZ_merged.bam" "${DIR}/bam/NPC_merged.bam" --labels SCZ NPC -p "${THREADS}" --plotFile "${DIR}/qc/fingerprint.png" --outRawCounts "${DIR}/qc/fingerprint_counts.tab" 2>> "${DIR}/logs/fingerprint.log"
multiBamSummary bins -b "${DIR}/bam/"*filtered.bam -p "${THREADS}" -o "${DIR}/qc/bam_summary.npz" 2>> "${DIR}/logs/multibamsummary.log"
plotCorrelation -in "${DIR}/qc/bam_summary.npz" --corMethod pearson --skipZeros --whatToPlot heatmap --colorMap RdYlBu --plotNumbers -o "${DIR}/qc/correlation_heatmap.png" 2>> "${DIR}/logs/plotcorrelation.log"

echo "======================================"
echo " GSE215991 COMPLETE: $(date)"
echo " Peaks: ${DIR}/peaks/"
echo "======================================"
