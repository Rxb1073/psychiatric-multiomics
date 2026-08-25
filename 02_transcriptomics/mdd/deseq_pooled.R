#!/usr/bin/env Rscript
suppressMessages({ library(DESeq2); library(sva) })
setwd("/rds/projects/e/elhamsak-dissertation-autism-2/MDD/transcriptomics/DE")
col <- read.delim("rna_colData.tsv", check.names=FALSE)
rownames(col) <- col$`sample-id`
cnt <- read.delim("merged_counts.tsv", row.names=1, check.names=FALSE)
cnt <- round(as.matrix(cnt[, rownames(col)]))
stopifnot(all(colnames(cnt)==rownames(col)))
col$disease <- factor(col$disease, levels=c("Control","MDD"))
col$study   <- factor(col$study)
cat("pooled:", ncol(cnt), "samples,", nrow(cnt), "genes\n")
print(table(col$study, col$disease))
keep <- rowSums(cnt >= 10) >= 0.5*ncol(cnt)
cnt <- cnt[keep, ]
cat("genes after filter:", nrow(cnt), "\n")
## ROUTE 1: ComBat-seq
adj <- ComBat_seq(cnt, batch=col$study, group=col$disease)
dds1 <- DESeqDataSetFromMatrix(adj, col, design=~ disease)
dds1 <- DESeq(dds1)
r1 <- results(dds1, contrast=c("disease","MDD","Control"), alpha=0.05)
r1 <- lfcShrink(dds1, coef="disease_MDD_vs_Control", res=r1, type="apeglm")
r1 <- as.data.frame(r1); r1$gene <- rownames(r1)
r1 <- r1[order(r1$padj), c("gene","log2FoldChange","lfcSE","pvalue","padj")]
write.table(r1, "deseq_pooled_ComBat_MDDvsControl.tsv", sep="\t", quote=FALSE, row.names=FALSE)
cat("\n=== POOLED (ComBat-seq): significant padj<0.05:", sum(r1$padj<0.05, na.rm=TRUE), "===\n")
## ROUTE 2: covariate sensitivity
dds2 <- DESeqDataSetFromMatrix(cnt, col, design=~ study + disease)
dds2 <- DESeq(dds2)
r2 <- results(dds2, contrast=c("disease","MDD","Control"), alpha=0.05)
r2 <- lfcShrink(dds2, coef="disease_MDD_vs_Control", res=r2, type="apeglm")
r2 <- as.data.frame(r2); r2$gene <- rownames(r2)
r2 <- r2[order(r2$padj), c("gene","log2FoldChange","lfcSE","pvalue","padj")]
write.table(r2, "deseq_pooled_covariate_MDDvsControl.tsv", sep="\t", quote=FALSE, row.names=FALSE)
cat("=== POOLED (covariate ~study+disease): significant padj<0.05:", sum(r2$padj<0.05, na.rm=TRUE), "===\n")
s1 <- r1$gene[which(r1$padj<0.05)]; s2 <- r2$gene[which(r2$padj<0.05)]
cat("\nComBat sig:", length(s1), " covariate sig:", length(s2), " overlap:", length(intersect(s1,s2)), "\n")
saveRDS(list(combat=r1, covariate=r2), "deseq_pooled_results.rds")
