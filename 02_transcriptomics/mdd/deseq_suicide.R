#!/usr/bin/env Rscript
suppressMessages({ library(DESeq2) })
setwd("/rds/projects/e/elhamsak-dissertation-autism-2/MDD/transcriptomics/DE")

## GSE101521 only — the three-level design (Control / MDD_suicide / MDD_nonsuicide)
col <- read.delim("GSE101521_colData.tsv", check.names=FALSE)
rownames(col) <- col$`sample-id`
col <- col[col$group3 != "NA", ]                 # keep the labelled samples
cnt <- read.delim("merged_counts.tsv", row.names=1, check.names=FALSE)
cnt <- round(as.matrix(cnt[, rownames(col)]))
stopifnot(all(colnames(cnt)==rownames(col)))

col$group3 <- factor(col$group3, levels=c("Control","MDD_nonsuicide","MDD_suicide"))
cat("group3 counts:\n"); print(table(col$group3))

## covariates present+complete for GSE101521 were sex, ph (from per-dataset run)
usecov <- c()
for (cv in c("sex","ph")) {
  v <- col[[cv]]
  if (all(!is.na(v)) && all(v!="") && length(unique(v))>1) {
    num <- suppressWarnings(as.numeric(v))
    col[[cv]] <- if (all(!is.na(num))) num else factor(v)
    usecov <- c(usecov, cv)
  }
}
form <- as.formula(paste("~", paste(c(usecov,"group3"), collapse=" + ")))
cat("design:", deparse(form), "\n")

dds <- DESeqDataSetFromMatrix(cnt, col, design=form)
dds <- dds[rowSums(counts(dds) >= 10) >= 0.5*ncol(dds), ]
cat("genes after filter:", nrow(dds), "\n")
dds <- DESeq(dds)

## three contrasts
contrasts <- list(
  MDDsuicide_vs_Control    = c("group3","MDD_suicide","Control"),
  MDDnonsuicide_vs_Control = c("group3","MDD_nonsuicide","Control"),
  suicide_vs_nonsuicide    = c("group3","MDD_suicide","MDD_nonsuicide")
)
for (nm in names(contrasts)) {
  r <- results(dds, contrast=contrasts[[nm]], alpha=0.05)
  r <- as.data.frame(r); r$gene <- rownames(r)
  r <- r[order(r$padj), c("gene","log2FoldChange","lfcSE","pvalue","padj")]
  write.table(r, paste0("deseq_suicide_",nm,".tsv"), sep="\t", quote=FALSE, row.names=FALSE)
  nsig <- sum(r$padj<0.05 & abs(r$log2FoldChange)>0.3, na.rm=TRUE)
  cat(sprintf("\n=== %s: sig (padj<0.05 & |lfc|>0.3): %d ===\n", nm, nsig))
  top <- r[!is.na(r$padj) & r$padj<0.05 & abs(r$log2FoldChange)>0.3, ]
  if (nrow(top)) print(head(top,10), row.names=FALSE)
}
