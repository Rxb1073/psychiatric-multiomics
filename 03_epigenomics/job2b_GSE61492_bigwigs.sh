#!/bin/bash
#SBATCH --job-name=GSE61492_BW
#SBATCH --account=elhamsak-dissertation-autism
#SBATCH --qos=bbdefault
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=/rds/projects/e/elhamsak-dissertation-autism/Epigenomics/logs/GSE61492_BW_%j.out
#SBATCH --error=/rds/projects/e/elhamsak-dissertation-autism/Epigenomics/logs/GSE61492_BW_%j.err

set -euo pipefail
THREADS=${SLURM_CPUS_PER_TASK}

module purge
module load bear-apps/2022a/live
module load MACS2/2.2.9.1-foss-2022a
module load deepTools/3.5.2-foss-2022a

DIR="/rds/projects/e/elhamsak-dissertation-autism/Epigenomics/epigenomics-CHD8/GSE61492"
INPUT_SAMPLE="SRR1576805"
CHIP_SAMPLES=(SRR1576798 SRR1576799 SRR1576800 SRR1576801 SRR1576802 SRR1576803 SRR1576804)

echo "[GSE61492] BigWigs..."
for SRR in "${CHIP_SAMPLES[@]}"; do
    bamCompare \
        -b1 "${DIR}/bam/${SRR}.filtered.bam" \
        -b2 "${DIR}/bam/${INPUT_SAMPLE}.filtered.bam" \
        --binSize 10 --normalizeUsing RPKM \
        --operation log2 --scaleFactorsMethod None \
        -p "${THREADS}" \
        -o "${DIR}/bigwig/${SRR}_log2FC.bw" \
        2>> "${DIR}/logs/${SRR}_bamcompare.log"
    echo "  Done: ${SRR}"
done

echo "[GSE61492] QC..."
BAMS=()
for SRR in "${CHIP_SAMPLES[@]}" "${INPUT_SAMPLE}"; do
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
echo " GSE61492 BigWigs + QC COMPLETE: $(date)"
echo "======================================"
