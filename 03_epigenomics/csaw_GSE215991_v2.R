
suppressPackageStartupMessages({library(csaw); library(edgeR); library(rtracklayer)})
DIR <- "/rds/projects/e/elhamsak-dissertation-autism/Epigenomics/epigenomics-ChIP-seq/GSE215991"
OUT <- file.path(DIR, "csaw")

data   <- readRDS(file.path(OUT, "windows_raw.rds"))
binned <- readRDS(file.path(OUT, "bins_10kb.rds"))
grp    <- factor(c(rep("SCZ",15), rep("NPC",15)), levels=c("NPC","SCZ"))

## --- diagnostics: why did the old filter fail? ---
fstat <- filterWindowsGlobal(data, binned)
message("enrichment-over-background quantiles:")
print(round(quantile(fstat$filter, c(.5,.9,.95,.99,.999,1)), 3))
message("  n > log2(3): ", sum(fstat$filter > log2(3)))
message("  n > log2(1.5): ", sum(fstat$filter > log2(1.5)))
message("  n > 0: ", sum(fstat$filter > 0))

## --- abundance filter (correct for broad marks) ---
ab <- aveLogCPM(asDGEList(data))
message("\nabundance (aveLogCPM) quantiles:")
print(round(quantile(ab, c(.1,.25,.5,.75,.9)), 3))
thr  <- aveLogCPM(10, lib.size=mean(data$totals))   # >= ~10 reads/window
keep <- ab > thr
message("\nkeeping windows with >= ~10 reads: ", sum(keep), " / ", length(keep),
        " (", round(100*mean(keep),1), "%)")
filtered <- data[keep, ]

filtered <- normFactors(binned, se.out=filtered)
y <- asDGEList(filtered)
design <- model.matrix(~0 + grp); colnames(design) <- levels(grp)
y   <- estimateDisp(y, design)
fit <- glmQLFit(y, design, robust=TRUE)
res <- glmQLFTest(fit, contrast=makeContrasts(SCZ - NPC, levels=design))
message("common BCV: ", round(sqrt(y$common.dispersion), 4))

merged  <- mergeWindows(rowRanges(filtered), tol=5000, max.width=100000)
tab     <- combineTests(merged$ids, res$table)
regions <- merged$regions; mcols(regions) <- tab
sig     <- regions[tab$FDR < 0.05]
message("regions tested: ", length(regions), " | FDR<0.05: ", length(sig))
if (length(sig)) {
  message("  up in SCZ: ", sum(sig$direction=="up"), " | down: ", sum(sig$direction=="down"))
  message("  median width: ", round(median(width(sig))), " bp")
}

write.csv(as.data.frame(regions), file.path(OUT,"H3K9me2_regions_abundance_filter.csv"), row.names=FALSE)
if (length(sig)) {
  write.csv(as.data.frame(sig), file.path(OUT,"H3K9me2_SCZ_vs_NPC_FDR05_v2.csv"), row.names=FALSE)
  export(sig, file.path(OUT,"H3K9me2_diff_regions_v2.bed"), format="BED")
}
save(filtered, y, fit, res, merged, tab, regions, file=file.path(OUT,"GSE215991_csaw_v2.RData"))

pdf(file.path(OUT,"csaw_qc_v2.pdf"), width=11, height=8)
par(mfrow=c(2,2))
hist(fstat$filter, breaks=100, main="log2 enrichment over background\n(near zero = broad mark)", xlab="log2 FC")
abline(v=log2(3), col="red", lwd=2)
hist(ab, breaks=100, main="abundance (aveLogCPM)", xlab="aveLogCPM"); abline(v=thr, col="red", lwd=2)
plotBCV(y)
plotMDS(cpm(y, log=TRUE), col=c("steelblue","tomato")[grp], pch=16, main="MDS (red=SCZ, blue=NPC)")
dev.off()
message("=== v2 COMPLETE ===")
