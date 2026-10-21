
suppressPackageStartupMessages({
  library(csaw); library(edgeR); library(GenomicRanges); library(rtracklayer)
})

DIR <- "/rds/projects/e/elhamsak-dissertation-autism/Epigenomics/epigenomics-ChIP-seq/GSE215991"
OUT <- file.path(DIR, "csaw"); dir.create(OUT, showWarnings=FALSE, recursive=TRUE)

SCZ <- sprintf("SRR219486%02d", 48:62)
NPC <- sprintf("SRR219486%02d", 63:77)
bam <- file.path(DIR, "bam", paste0(c(SCZ, NPC), ".filtered.bam"))
stopifnot(all(file.exists(bam)), all(file.exists(paste0(bam, ".bai"))))
grp <- factor(c(rep("SCZ", 15), rep("NPC", 15)), levels=c("NPC","SCZ"))
message("BAMs: ", length(bam), " | SCZ=", sum(grp=="SCZ"), " NPC=", sum(grp=="NPC"))

# BAMs are already MAPQ>=30 + blacklist filtered; restrict to main chromosomes
param <- readParam(pe="none", restrict=paste0("chr", c(1:22,"X","Y")))

## ---- 1. Fragment length (cross-correlation on chr1, 3 samples) ----
message("[1] estimating fragment length ...")
p1 <- reform(param, restrict="chr1")
fl <- sapply(bam[c(1, 8, 16)], function(b) maximizeCcf(correlateReads(b, max.dist=500, param=p1)))
EXT <- round(median(fl)); message("    per-sample: ", paste(fl, collapse=", "), " | using ext=", EXT)

## ---- 2. Count reads in 2 kb windows ----
message("[2] windowCounts 2kb (slow) ...")
data <- windowCounts(bam, width=2000, spacing=2000, ext=EXT, param=param)
colData(data)$group <- grp
message("    windows: ", nrow(data))
saveRDS(data, file.path(OUT, "windows_raw.rds"))

## ---- 3. Large bins: background model + normalisation basis ----
message("[3] binning 10kb ...")
binned <- windowCounts(bam, bin=TRUE, width=10000, param=param)
saveRDS(binned, file.path(OUT, "bins_10kb.rds"))

## ---- 4. Filter windows vs global background ----
message("[4] filtering vs global background ...")
fstat <- filterWindowsGlobal(data, binned)
keep  <- fstat$filter > log2(3)          # >=3x above background
message("    kept ", sum(keep), " / ", length(keep), " (", round(100*mean(keep),1), "%)")
filtered <- data[keep, ]

## ---- 5. TMM normalisation from the 10kb bins (composition bias) ----
message("[5] TMM on bins -> applied to windows ...")
filtered <- normFactors(binned, se.out=filtered)
print(data.frame(sample=c(SCZ,NPC), group=grp,
                 norm.factor=round(filtered$norm.factors,3)))

## ---- 6. edgeR quasi-likelihood, SCZ vs NPC ----
message("[6] edgeR QL ...")
y <- asDGEList(filtered)
design <- model.matrix(~0 + grp); colnames(design) <- levels(grp)
y   <- estimateDisp(y, design)
fit <- glmQLFit(y, design, robust=TRUE)
res <- glmQLFTest(fit, contrast=makeContrasts(SCZ - NPC, levels=design))
message("    common BCV: ", round(sqrt(y$common.dispersion), 4))

## ---- 7. Merge adjacent windows into regions ----
message("[7] merging windows ...")
merged <- mergeWindows(rowRanges(filtered), tol=1000, max.width=10000)
tab    <- combineTests(merged$ids, res$table)
regions <- merged$regions
mcols(regions) <- tab

sig <- regions[tab$FDR < 0.05]
message("    regions: ", length(regions), " | FDR<0.05: ", length(sig))
if (length(sig)) {
  message("    up in SCZ: ",   sum(sig$direction == "up"))
  message("    down in SCZ: ", sum(sig$direction == "down"))
}

## ---- 8. Save ----
df <- as.data.frame(regions)
write.csv(df, file.path(OUT, "H3K9me2_all_regions.csv"), row.names=FALSE)
write.csv(df[df$FDR < 0.05, ], file.path(OUT, "H3K9me2_SCZ_vs_NPC_FDR05.csv"), row.names=FALSE)
if (length(sig)) export(sig, file.path(OUT, "H3K9me2_diff_regions.bed"), format="BED")
save(filtered, y, fit, res, merged, tab, regions, EXT,
     file=file.path(OUT, "GSE215991_csaw.RData"))

pdf(file.path(OUT, "csaw_qc.pdf"), width=11, height=8)
par(mfrow=c(2,2))
hist(fstat$filter, breaks=100, main="log2 enrichment over background", xlab="log2 FC")
abline(v=log2(3), col="red", lwd=2)
plotBCV(y)
plotMDS(cpm(y, log=TRUE), col=c("tomato","steelblue")[grp], pch=16,
        main="MDS - H3K9me2 (red=SCZ, blue=NPC)")
plotQLDisp(fit)
dev.off()

message("=== csaw COMPLETE: ", length(sig), " differential regions ===")
