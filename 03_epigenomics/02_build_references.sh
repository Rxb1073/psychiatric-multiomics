#!/bin/bash
#SBATCH --job-name=build_references
#SBATCH --account=elhamsak-dissertation-autism
#SBATCH --qos=bbdefault
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=/rds/projects/e/elhamsak-dissertation-autism/Epigenomics/logs/build_refs_%j.out
#SBATCH --error=/rds/projects/e/elhamsak-dissertation-autism/Epigenomics/logs/build_refs_%j.err
 
# =============================================================================
#  02_build_references.sh
#  Downloads blacklists + builds Bowtie2 indices for hg19 and hg38
#  Run ONCE before submitting the 6 analysis jobs
# =============================================================================
 
set -euo pipefail
THREADS=${SLURM_CPUS_PER_TASK}
 
REFS="/rds/projects/e/elhamsak-dissertation-autism/Epigenomics/refs"
MIN_SPACE_GB=20
 
echo "=============================================="
echo "  BUILD REFERENCES"
echo "  Job:   ${SLURM_JOB_ID}"
echo "  Start: $(date)"
echo "=============================================="
 
# ── Load modules (2023a for Bowtie2) ──────────────────────────────────────
module purge
module load bear-apps/2023a/live
module load Bowtie2/2.5.4-GCC-12.3.0
module load SAMtools/1.21-GCC-12.3.0
 
# ── Space check ────────────────────────────────────────────────────────────
AVAIL_GB=$(df -BG "${REFS%/*}" | awk 'NR==2 {gsub("G","",$4); print $4}')
echo "  Available: ${AVAIL_GB}GB | Required: ${MIN_SPACE_GB}GB"
if [[ "${AVAIL_GB}" -lt "${MIN_SPACE_GB}" ]]; then
    echo "  ✗ NOT ENOUGH SPACE — stopping"
    exit 1
fi
echo "  ✓ Space check passed"
 
# ── STEP 1: Create folders ─────────────────────────────────────────────────
echo ""
echo "STEP 1: Creating folders..."
mkdir -p "${REFS}"/{logs,bowtie2/hg19,bowtie2/hg38,blacklists,genomes/hg19,genomes/hg38}
echo "  ✓ Folders created"
 
# ── STEP 2: Blacklists ────────────────────────────────────────────────────
echo ""
echo "STEP 2: Downloading blacklists..."
 
if [[ ! -f "${REFS}/blacklists/hg19-blacklist.v2.bed" ]]; then
    wget -q --show-progress \
        -O "${REFS}/blacklists/hg19-blacklist.v2.bed.gz" \
        "https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg19-blacklist.v2.bed.gz"
    gunzip "${REFS}/blacklists/hg19-blacklist.v2.bed.gz"
    echo "  ✓ hg19 blacklist"
else
    echo "  Skip hg19 blacklist (exists)"
fi
 
if [[ ! -f "${REFS}/blacklists/hg38-blacklist.v2.bed" ]]; then
    wget -q --show-progress \
        -O "${REFS}/blacklists/hg38-blacklist.v2.bed.gz" \
        "https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz"
    gunzip "${REFS}/blacklists/hg38-blacklist.v2.bed.gz"
    echo "  ✓ hg38 blacklist"
else
    echo "  Skip hg38 blacklist (exists)"
fi
 
# ── STEP 3: hg19 genome + index ───────────────────────────────────────────
echo ""
echo "STEP 3: hg19 genome + Bowtie2 index..."
HG19_FA="${REFS}/genomes/hg19/hg19.fa"
HG19_IDX="${REFS}/bowtie2/hg19/hg19"
 
if [[ ! -f "${HG19_IDX}.1.bt2" ]]; then
    if [[ ! -f "${HG19_FA}" ]]; then
        echo "  Downloading hg19 genome (~3GB)..."
        wget -q --show-progress \
            -O "${HG19_FA}.gz" \
            "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz"
        gunzip "${HG19_FA}.gz"
        echo "  ✓ hg19 genome downloaded"
    fi
    echo "  Building hg19 index (~1 hour)..."
    bowtie2-build --threads "${THREADS}" "${HG19_FA}" "${HG19_IDX}" \
        2>&1 | tee "${REFS}/logs/bowtie2_build_hg19.log"
    echo "  ✓ hg19 index built"
else
    echo "  Skip hg19 index (exists)"
fi
 
# ── STEP 4: hg38 genome + index ───────────────────────────────────────────
echo ""
echo "STEP 4: hg38 genome + Bowtie2 index..."
HG38_FA="${REFS}/genomes/hg38/hg38.fa"
HG38_IDX="${REFS}/bowtie2/hg38/hg38"
 
if [[ ! -f "${HG38_IDX}.1.bt2" ]]; then
    if [[ ! -f "${HG38_FA}" ]]; then
        echo "  Downloading hg38 genome (~3.5GB)..."
        wget -q --show-progress \
            -O "${HG38_FA}.gz" \
            "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"
        gunzip "${HG38_FA}.gz"
        echo "  ✓ hg38 genome downloaded"
    fi
    echo "  Building hg38 index (~1 hour)..."
    bowtie2-build --threads "${THREADS}" "${HG38_FA}" "${HG38_IDX}" \
        2>&1 | tee "${REFS}/logs/bowtie2_build_hg38.log"
    echo "  ✓ hg38 index built"
else
    echo "  Skip hg38 index (exists)"
fi
 
# ── STEP 5: Install R packages ────────────────────────────────────────────
echo ""
echo "STEP 5: Installing R packages..."
module purge
module load bear-apps/2023a/live
module load R/4.4.1-gfbf-2023a
 
Rscript - << 'RSCRIPT'
message("Installing R packages...")
lib_path <- Sys.getenv("R_LIBS_USER")
if (lib_path == "") lib_path <- "~/R/library"
dir.create(lib_path, recursive=TRUE, showWarnings=FALSE)
.libPaths(lib_path)
 
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager", repos="https://cloud.r-project.org", lib=lib_path)
library(BiocManager)
 
pkgs <- c(
    "minfi", "limma", "DMRcate", "missMethyl", "bumphunter",
    "IlluminaHumanMethylation450kanno.ilmn12.hg19",
    "IlluminaHumanMethylation450kmanifest",
    "GEOquery", "ggplot2", "pheatmap", "sva"
)
installed  <- rownames(installed.packages(lib=lib_path))
to_install <- pkgs[!pkgs %in% installed]
if (length(to_install) == 0) {
    message("All packages already installed!")
} else {
    message("Installing ", length(to_install), " packages...")
    BiocManager::install(to_install, lib=lib_path, ask=FALSE, update=FALSE)
}
failed <- c()
for (pkg in pkgs) {
    ok <- requireNamespace(pkg, quietly=TRUE)
    message("  ", ifelse(ok,"✓","✗"), " ", pkg)
    if (!ok) failed <- c(failed, pkg)
}
if (length(failed) > 0) {
    message("WARNING: Failed: ", paste(failed, collapse=", "))
    quit(status=1)
}
message("All R packages installed!")
RSCRIPT
 
# ── STEP 6: Verify ────────────────────────────────────────────────────────
echo ""
echo "STEP 6: Verification..."
ERRORS=0
 
for BED in \
    "${REFS}/blacklists/hg19-blacklist.v2.bed" \
    "${REFS}/blacklists/hg38-blacklist.v2.bed"; do
    [[ -f "${BED}" ]] && echo "  ✓ $(basename ${BED})" || { echo "  ✗ MISSING: ${BED}"; ERRORS=$((ERRORS+1)); }
done
 
for IDX in "${REFS}/bowtie2/hg19/hg19" "${REFS}/bowtie2/hg38/hg38"; do
    BT2_COUNT=$(ls "${IDX}".*.bt2 2>/dev/null | wc -l)
    if [[ "${BT2_COUNT}" -eq 6 ]]; then
        echo "  ✓ Bowtie2 index: $(basename ${IDX})"
    else
        echo "  ✗ INCOMPLETE: ${IDX} (${BT2_COUNT}/6 files)"
        ERRORS=$((ERRORS+1))
    fi
done
 
echo ""
if [[ "${ERRORS}" -eq 0 ]]; then
    echo "  ✓ ALL REFERENCES READY"
    echo ""
    echo "  Reference paths:"
    echo "  BOWTIE2_HG19=\"${REFS}/bowtie2/hg19/hg19\""
    echo "  BOWTIE2_HG38=\"${REFS}/bowtie2/hg38/hg38\""
    echo "  BLACKLIST_HG19=\"${REFS}/blacklists/hg19-blacklist.v2.bed\""
    echo "  BLACKLIST_HG38=\"${REFS}/blacklists/hg38-blacklist.v2.bed\""
else
    echo "  ✗ ${ERRORS} ERROR(S) — fix before running analysis jobs"
    exit 1
fi
 
echo ""
echo "=============================================="
echo "  COMPLETE: $(date)"
echo "  Next: bash submit_all_jobs.sh"
echo "=============================================="
 