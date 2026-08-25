#!/usr/bin/env Rscript
suppressMessages({ library(DESeq2) })
setwd("/rds/projects/e/elhamsak-dissertation-autism-2/MDD/transcriptomics/DE")
args <- commandArgs(trailingOnly=TRUE)
DSNAME <- args[1]
col <- read.delim(paste0(DSNAME,"_colData.tsv"), check.names=FALSE)
rownames(col) <- col$`sample-id`
cnt <- read.delim("merged_counts.tsv", row.names=1, check.names=FALSE)
cnt <- cnt[, rownames(col)]
stopifnot(all(colnames(cnt)==rownames(col)))
cat(DSNAME, ":", ncol(cnt), "samples,", nrow(cnt), "genes\n")
col$disease <- factor(col$disease, levels=c("Control","MDD"))
cand <- c("sex","age","pmi","rin","ph"); usecov <- c()
for (cv in cand) {
  v <- col[[cv]]
  if (!is.null(v) && all(!is.na(v)) && all(v!="") && length(unique(v))>1) {
    num <- suppressWarnings(as.numeric(v))
    if (all(!is.na(num))) { col[[cv]] <- num } else { col[[cv]] <- factor(v) }
    usecov <- c(usecov, cv)
  }
}
cat("covariates used:", if(length(usecov)) paste(usecov,collapse=", ") else "(none)", "\n")
form <- as.formula(paste("~", paste(c(usecov,"disease"), collapse=" + ")))
cat("design:", deparse(form), "\n")
dds <- DESeqDataSetFromMatrix(round(as.matrix(cnt)), col, design=form)
dds <- dds[rowSums(counts(dds) >= 10) >= 0.5*ncol(dds), ]
cat("genes after filter:", nrow(dds), "\n")
dds <- DESeq(dds)
res <- results(dds, contrast=c("disease","MDD","Control"), alpha=0.05)
res <- lfcShrink(dds, coef="disease_MDD_vs_Control", res=res, type="apeglm")
res <- as.data.frame(res); res$gene <- rownames(res)
res <- res[order(res$padj), c("gene","log2FoldChange","lfcSE","pvalue","padj")]
write.table(res, paste0("deseq_",DSNAME,"_MDDvsControl.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
nsig <- sum(res$padj < 0.05, na.rm=TRUE)
cat("\n===", DSNAME, ": significant genes (padj<0.05):", nsig, "===\n")
print(head(res[!is.na(res$padj) & res$padj<0.05, ], 15), row.names=FALSE)
