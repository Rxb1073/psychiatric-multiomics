.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths())); suppressMessages(library(DESeq2))
logf <- function(...) cat(sprintf(...),"\n")
NP <- "/rds/projects/e/elhamsak-dissertation-neuropsychiatric/RNA_seq"
SHD <- file.path(NP,"counts_scz/sheets")
OR <- "/rds/projects/e/elhamsak-dissertation-autism-2/RNA_seq/DE_analysis"
bt <- read.delim(file.path(OR,"gene_id2biotype.tsv"), header=FALSE, col.names=c("gid","bt"))
keep_bt <- bt$gid[bt$bt %in% c("protein_coding","lncRNA")]
cm <- read.delim(file.path(NP,"PRJNA574470/counts/PRJNA574470_gene_counts_matrix.tsv"),
                 comment.char="#", check.names=FALSE)
g <- cm[[1]]; cm[[1]] <- NULL
cm <- cm[, !colnames(cm) %in% c("Chr","Start","End","Strand","Length"), drop=FALSE]
colnames(cm) <- regmatches(colnames(cm), regexpr("SRR[0-9]+", colnames(cm)))
cnt <- as.matrix(cm); rownames(cnt) <- g; storage.mode(cnt) <- "integer"
ss <- merge(read.delim(file.path(SHD,"GSE138082_samplesheet.tsv"), check.names=FALSE),
            read.delim(file.path(SHD,"GSE138082_prua.tsv")), by="srr")
rownames(ss) <- ss$srr
k <- intersect(colnames(cnt), rownames(ss)); cnt <- cnt[,k,drop=FALSE]; ss <- ss[k,]
cnt <- cnt[rownames(cnt) %in% keep_bt,,drop=FALSE]
run <- function(ss,cnt,form,lab){
  ss$diagnosis <- relevel(factor(ss$diagnosis),"Control")
  for (v in c("age","pmi","rin","prua")) ss[[paste0(v,"_z")]] <- as.numeric(scale(as.numeric(ss[[v]])))
  cc <- cnt[rowSums(cnt>=10) >= min(table(ss$diagnosis)),,drop=FALSE]
  d <- suppressMessages(DESeq(DESeqDataSetFromMatrix(cc,ss,form), quiet=TRUE))
  r <- as.data.frame(results(d, name=tail(resultsNames(d),1)))
  s <- function(l) sum(!is.na(r$padj)&r$padj<0.05&abs(r$log2FoldChange)>=l)
  logf("  %-34s n=%d  padj05:%5d | >=0.3:%5d | >=1.0:%4d", lab, nrow(ss), s(0), s(0.3), s(1))
  invisible(rownames(r)[!is.na(r$padj)&r$padj<0.05])
}
logf("=== ALL 78 ===")
a1 <- run(ss,cnt, ~age_z+sex+pmi_z+rin_z+diagnosis,          "R1 +RIN")
a2 <- run(ss,cnt, ~age_z+sex+pmi_z+rin_z+prua_z+diagnosis,   "R1b +RIN +PRUA")
logf("  R1 -> R1b retained: %d/%d", length(intersect(a1,a2)), length(a1))
logf("\n=== EXCLUDING SRR10192704 (PRUA 25.3%%) ===")
kk <- rownames(ss)!="SRR10192704"
b1 <- run(ss[kk,],cnt[,kk], ~age_z+sex+pmi_z+rin_z+diagnosis,        "R1 +RIN")
b2 <- run(ss[kk,],cnt[,kk], ~age_z+sex+pmi_z+rin_z+prua_z+diagnosis, "R1b +RIN +PRUA")
logf("  R1 full -> R1 n77 retained: %d/%d", length(intersect(a1,b1)), length(a1))
logf("\n=== EXCLUDING PRUA < 65%% (5 samples) ===")
k5 <- ss$prua >= 65
c1 <- run(ss[k5,],cnt[,k5], ~age_z+sex+pmi_z+rin_z+diagnosis, "R1 +RIN")
logf("  R1 full -> R1 n%d retained: %d/%d", sum(k5), length(intersect(a1,c1)), length(a1))
