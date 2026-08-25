#!/usr/bin/env Rscript
# DESeq2 per-dataset DE — SCZ vs Control. Adapted from ASD v3.
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
suppressMessages(library(DESeq2))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Provide a dataset key")
DS <- args[1]; CASE <- "SCZ"
logf <- function(...) cat(sprintf(...), "\n")
NEURO   <- "/rds/projects/e/elhamsak-dissertation-neuropsychiatric/RNA_seq"
CS      <- file.path(NEURO, "counts_scz")
SHEETS  <- file.path(CS, "sheets")
OUTROOT <- file.path(CS, "DE_per_dataset")
ASD_DE  <- "/rds/projects/e/elhamsak-dissertation-autism-2/RNA_seq/DE_analysis"
GTFMAP  <- file.path(ASD_DE, "gene_id2name.tsv")
BIOTYPE <- file.path(ASD_DE, "gene_id2biotype.tsv")
KEEP_BIOTYPES <- c("protein_coding", "lncRNA")
dir.create(OUTROOT, recursive = TRUE, showWarnings = FALSE)
cfg <- list(
  GSE80655_DLPFC = list(
    counts=file.path(NEURO,"PRJNA319583/counts/PRJNA319583_gene_counts_matrix.tsv"),
    sheet=file.path(SHEETS,"GSE80655_DLPFC_samplesheet.tsv"),
    design="~ age + sex + pmi + ph + diagnosis",
    numeric_vars=c("age","pmi","ph"),
    note="DLPFC; full covariates; 28 SCZ / 29 Control"),
  GSE78936 = list(
    counts=file.path(NEURO,"PRJNA314463/counts/PRJNA314463_gene_counts.txt"),
    sheet=file.path(SHEETS,"GSE78936_samplesheet.tsv"),
    design="~ brain_region + diagnosis",
    numeric_vars=character(0),
    note="frontal BA9/11/24; subregion only; 28 SCZ / 24 Control"),
  GSE138082 = list(
    counts=file.path(NEURO,"PRJNA574470/counts/PRJNA574470_gene_counts_matrix.tsv"),
    sheet=file.path(SHEETS,"GSE138082_samplesheet.tsv"),
    design="~ age + sex + pmi + rin + diagnosis",
    numeric_vars=c("age","pmi","rin"),
    note="hippocampus; full covariates; 39 SCZ / 39 Control")
)
if (!DS %in% names(cfg)) stop(sprintf("Unknown dataset '%s'", DS))
C <- cfg[[DS]]; OUT <- file.path(OUTROOT, DS)
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
logf("=== DESeq2 [SCZ]: %s ===", DS); logf("note: %s", C$note); logf("design: %s", C$design)
logf("counts: %s", C$counts); logf("sheet : %s", C$sheet)
ss <- read.delim(C$sheet, stringsAsFactors=FALSE, check.names=FALSE)
stopifnot("srr" %in% names(ss), "diagnosis" %in% names(ss)); rownames(ss) <- ss$srr
cm <- read.delim(C$counts, comment.char="#", check.names=FALSE, stringsAsFactors=FALSE)
gid <- cm[[1]]; cm[[1]] <- NULL
cm <- cm[, !(colnames(cm) %in% c("Chr","Start","End","Strand","Length")), drop=FALSE]
srr_hit <- regmatches(colnames(cm), regexpr("SRR[0-9]+", colnames(cm)))
has_srr <- nchar(srr_hit) > 0
if (!all(has_srr)) logf("WARN: %d cols no SRR, dropped", sum(!has_srr))
cm <- cm[, has_srr, drop=FALSE]
colnames(cm) <- regmatches(colnames(cm), regexpr("SRR[0-9]+", colnames(cm)))
rownames(cm) <- gid
cmat <- as.matrix(cm); cmat[is.na(cmat)] <- 0; storage.mode(cmat) <- "integer"
logf("counts matrix: %d genes x %d samples", nrow(cmat), ncol(cmat))
if (file.exists(BIOTYPE)) {
  bt <- read.delim(BIOTYPE, header=FALSE, col.names=c("gid","biotype"), stringsAsFactors=FALSE)
  btv <- setNames(bt$biotype, bt$gid); before <- nrow(cmat); gb <- btv[rownames(cmat)]
  logf("--- biotype top6 ---"); print(head(sort(table(gb, useNA="ifany"), decreasing=TRUE),6))
  cmat <- cmat[!is.na(gb) & gb %in% KEEP_BIOTYPES, , drop=FALSE]
  logf("biotype filter: %d -> %d genes", before, nrow(cmat))
} else logf("WARN: no biotype file!")
common <- intersect(colnames(cmat), rownames(ss))
logf("samples in both: %d (sheet %d, counts %d)", length(common), nrow(ss), ncol(cmat))
miss <- setdiff(rownames(ss), colnames(cmat))
if (length(miss)) logf("WARN: %d sheet samples missing: %s", length(miss), paste(miss, collapse=","))
cmat <- cmat[, common, drop=FALSE]; ss <- ss[common, , drop=FALSE]
stopifnot(identical(colnames(cmat), rownames(ss)))
design_vars <- all.vars(as.formula(C$design))
for (v in design_vars) {
  if (!v %in% names(ss)) stop(sprintf("design var '%s' not in sheet", v))
  if (v %in% C$numeric_vars) ss[[v]] <- as.numeric(ss[[v]]) else ss[[v]] <- factor(ss[[v]])
}
ss$diagnosis <- relevel(factor(ss$diagnosis), ref="Control")
na_rows <- Reduce(`|`, lapply(design_vars, function(v) is.na(ss[[v]])))
if (any(na_rows)) { logf("WARN: dropping %d NA-covariate samples", sum(na_rows))
  ss <- ss[!na_rows, , drop=FALSE]; cmat <- cmat[, rownames(ss), drop=FALSE] }
for (v in C$numeric_vars) {
  mu <- mean(ss[[v]]); sdv <- stats::sd(ss[[v]])
  logf("centering '%s': mean=%.3f sd=%.3f (raw %.2f-%.2f)", v, mu, sdv, min(ss[[v]]), max(ss[[v]]))
  ss[[v]] <- if (sdv>0) as.numeric((ss[[v]]-mu)/sdv) else as.numeric(ss[[v]]-mu)
}
logf("diagnosis table:"); print(table(ss$diagnosis))
for (v in setdiff(design_vars, c("diagnosis", C$numeric_vars))) {
  cat(sprintf("  %s x diagnosis:\n", v)); print(table(ss[[v]], ss$diagnosis)) }
dds <- DESeqDataSetFromMatrix(cmat, colData=ss, design=as.formula(C$design))
sg <- min(table(ss$diagnosis)); keep <- rowSums(counts(dds) >= 10) >= sg
logf("genes passing filter: %d / %d (smallest group n=%d)", sum(keep), nrow(dds), sg)
dds <- dds[keep, ]; dds <- DESeq(dds)
logf("resultsNames: %s", paste(resultsNames(dds), collapse=", "))
coef <- grep(sprintf("diagnosis_%s_vs_Control", CASE), resultsNames(dds), value=TRUE)
stopifnot(length(coef)==1)
res <- results(dds, contrast=c("diagnosis",CASE,"Control"), alpha=0.05)
shrink_type <- if (requireNamespace("apeglm", quietly=TRUE)) "apeglm" else "normal"
if (shrink_type=="normal") logf("WARN: apeglm missing, using normal")
logf("lfcShrink type: %s", shrink_type)
resL <- lfcShrink(dds, coef=coef, type=shrink_type)
gbio <- setNames(rep(NA_character_, nrow(res)), rownames(res))
if (file.exists(BIOTYPE)) { b <- read.delim(BIOTYPE, header=FALSE, col.names=c("gid","biotype")); gbio[] <- setNames(b$biotype,b$gid)[rownames(res)] }
gname <- setNames(rep(NA_character_, nrow(res)), rownames(res))
if (file.exists(GTFMAP)) { m <- read.delim(GTFMAP, header=FALSE, col.names=c("gid","name")); gname[] <- setNames(m$name,m$gid)[rownames(res)] }
out <- data.frame(gene_id=rownames(res), gene_name=gname[rownames(res)], biotype=gbio[rownames(res)],
  baseMean=res$baseMean, log2FC_MLE=res$log2FoldChange, log2FC_shrunk=resL$log2FoldChange,
  lfcSE_shrunk=resL$lfcSE, stat=res$stat, pvalue=res$pvalue, padj=res$padj,
  check.names=FALSE, row.names=NULL)
out <- out[order(out$padj), ]
write.table(out, file.path(OUT, sprintf("%s_SCZ_DESeq2_full_results.tsv", DS)), sep="\t", quote=FALSE, row.names=FALSE)
tiers <- list("padj05"=subset(out,!is.na(padj)&padj<0.05),
  "padj05_absLFC0.3"=subset(out,!is.na(padj)&padj<0.05&abs(log2FC_shrunk)>0.3),
  "padj05_absLFC0.5"=subset(out,!is.na(padj)&padj<0.05&abs(log2FC_shrunk)>0.5),
  "padj05_absLFC1"=subset(out,!is.na(padj)&padj<0.05&abs(log2FC_shrunk)>1))
logf("--- DE gene counts by threshold ---")
for (nm in names(tiers)) { t <- tiers[[nm]]
  logf("  %-18s : %5d (up=%d, down=%d)", nm, nrow(t), sum(t$log2FC_shrunk>0), sum(t$log2FC_shrunk<0))
  write.table(t, file.path(OUT, sprintf("%s_SCZ_sig_%s.tsv", DS, nm)), sep="\t", quote=FALSE, row.names=FALSE) }
nc <- counts(dds, normalized=TRUE)
write.table(data.frame(gene_id=rownames(nc), nc, check.names=FALSE),
  file.path(OUT, sprintf("%s_SCZ_normalized_counts.tsv", DS)), sep="\t", quote=FALSE, row.names=FALSE)
vsd <- tryCatch(vst(dds, blind=FALSE), error=function(e) varianceStabilizingTransformation(dds, blind=FALSE))
pdf(file.path(OUT, sprintf("%s_SCZ_QC.pdf", DS)), width=7, height=6)
pca <- plotPCA(vsd, intgroup="diagnosis", returnData=TRUE); pv <- round(100*attr(pca,"percentVar"))
plot(pca$PC1, pca$PC2, col=as.integer(pca$diagnosis)+1, pch=19,
  xlab=sprintf("PC1 (%d%%)",pv[1]), ylab=sprintf("PC2 (%d%%)",pv[2]), main=sprintf("%s PCA",DS))
legend("topright", legend=levels(pca$diagnosis), col=seq_along(levels(pca$diagnosis))+1, pch=19)
DESeq2::plotMA(resL, main=sprintf("%s MA",DS), ylim=c(-4,4))
plotDispEsts(dds, main=sprintf("%s disp",DS))
hist(res$pvalue, breaks=50, col="grey", main=sprintf("%s pval hist",DS), xlab="pvalue")
dev.off()
sink(file.path(OUT, sprintf("%s_SCZ_sessionInfo.txt", DS))); print(sessionInfo()); sink()
logf("outputs -> %s", OUT); logf("=== done: %s ===", DS)
