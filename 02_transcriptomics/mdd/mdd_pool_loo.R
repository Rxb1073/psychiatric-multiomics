.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
suppressMessages({library(DESeq2)})
logf <- function(...) cat(sprintf(...),"\n")
DE <- "/rds/projects/e/elhamsak-dissertation-autism-2/MDD/transcriptomics/DE"
NP <- "/rds/projects/e/elhamsak-dissertation-neuropsychiatric/RNA_seq"

cnt <- read.delim(file.path(DE,"merged_counts.tsv"), row.names=1, check.names=FALSE)
col <- read.delim(file.path(DE,"rna_colData.tsv"), check.names=FALSE)
rownames(col) <- col[["sample-id"]]
col <- col[colnames(cnt),]
logf("loaded: %d genes x %d samples", nrow(cnt), ncol(cnt))

rt <- read.csv(file.path(NP,"metadata/PRJNA319583_SraRunTable.csv"), check.names=FALSE)
bio <- setNames(rt$BioSample, rt$Run)
newid <- ifelse(colnames(cnt) %in% names(bio), bio[colnames(cnt)], colnames(cnt))
cnt2 <- t(rowsum(t(as.matrix(cnt)), group=newid))
keep <- !duplicated(newid); col2 <- col[keep,]; rownames(col2) <- newid[keep]
col2 <- col2[colnames(cnt2),]
logf("after collapse: %d samples (was %d)", ncol(cnt2), ncol(cnt))
print(table(col2$study, col2$disease))

run_de <- function(C, M, tag) {
  M$disease <- relevel(factor(M$disease), "Control")
  M$study <- factor(M$study)
  for (v in c("age","pmi","rin","ph")) {
    x <- suppressWarnings(as.numeric(as.character(M[[v]])))
    if (sum(!is.na(x)) >= 0.5*length(x)) { x[is.na(x)] <- median(x,na.rm=TRUE); M[[v]] <- as.numeric(scale(x)) }
    else M[[v]] <- NULL
  }
  cov <- intersect(c("sex","age","pmi","ph"), names(M))
  cov <- cov[sapply(cov, function(v) length(unique(M[[v]]))>1)]
  fm <- as.formula(paste("~ study +", paste(cov, collapse=" + "), "+ disease"))
  keep <- rowSums(C >= 10) >= min(table(M$disease))
  dds <- DESeq(DESeqDataSetFromMatrix(C[keep,], M, fm), quiet=TRUE)
  res <- results(dds, name="disease_MDD_vs_Control")
  n <- sum(res$padj < 0.05, na.rm=TRUE)
  logf("  %-24s %s", tag, deparse(fm))
  logf("  %-24s genes FDR<0.05: %d", tag, n)
  rownames(res)[which(res$padj < 0.05)]
}

logf("\n=========== POOLED, donor-correct ===========")
full <- run_de(cnt2, col2, "all 3 cohorts")
logf("\n=========== LEAVE-ONE-DATASET-OUT ===========")
for (s in unique(col2$study)) {
  k <- col2$study != s
  g <- run_de(cnt2[,k], col2[k,], paste("drop", s))
  logf("      retained of full: %d / %d\n", length(intersect(g, full)), length(full))
}
