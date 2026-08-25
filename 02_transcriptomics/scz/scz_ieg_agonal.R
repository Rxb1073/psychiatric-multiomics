## Are the IEG effects diagnosis, or agonal state (pH/PMI)?
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
suppressMessages({library(DESeq2); library(AnnotationDbi); library(org.Hs.eg.db)})
logf <- function(...) cat(sprintf(...),"\n")
NP <- "/rds/projects/e/elhamsak-dissertation-neuropsychiatric/RNA_seq"
CS <- file.path(NP,"counts_scz"); SHD <- file.path(CS,"sheets")
OR <- "/rds/projects/e/elhamsak-dissertation-autism-2/RNA_seq/DE_analysis"
RT <- "/rds/projects/e/elhamsak-dissertation-autism-2/MDD/transcriptomics/GSE80655/metadata/PRJNA319583_SraRunTable.csv"
IEG <- c("EGR1","EGR2","EGR3","NR4A1","NR4A2","NR4A3","FOS","FOSB","JUN","JUNB",
         "ARC","ATF3","DUSP1","DUSP6","NPAS4")

bt <- read.delim(file.path(OR,"gene_id2biotype.tsv"), header=FALSE, col.names=c("gid","bt"))
cm <- read.delim(file.path(NP,"PRJNA319583/counts/PRJNA319583_gene_counts_matrix.tsv"),
                 comment.char="#", check.names=FALSE)
g <- cm[[1]]; cm[[1]] <- NULL
cm <- cm[, !colnames(cm) %in% c("Chr","Start","End","Strand","Length"), drop=FALSE]
colnames(cm) <- regmatches(colnames(cm), regexpr("SRR[0-9]+", colnames(cm)))
cnt <- as.matrix(cm); rownames(cnt) <- g; storage.mode(cnt) <- "integer"
ss <- read.delim(file.path(SHD,"GSE80655_DLPFC_samplesheet.tsv"), check.names=FALSE)
ss <- merge(ss, read.delim(file.path(SHD,"GSE80655_DLPFC_prua.tsv"))[,c("srr","prua")], by="srr")
rownames(ss) <- ss$srr
k <- intersect(colnames(cnt), rownames(ss)); cnt <- cnt[,k]; ss <- ss[k,]
rt <- read.csv(RT, check.names=FALSE); bio <- setNames(rt$BioSample, rt$Run)[colnames(cnt)]
idx <- split(seq_len(ncol(cnt)), bio)
c2 <- vapply(idx, function(i) as.integer(rowSums(cnt[,i,drop=FALSE])), integer(nrow(cnt)))
rownames(c2) <- rownames(cnt); ke <- !duplicated(bio)
ss <- ss[ke,]; rownames(ss) <- bio[ke]; cnt <- c2[, rownames(ss)]
cnt <- cnt[rownames(cnt) %in% bt$gid[bt$bt %in% c("protein_coding","lncRNA")],]
ss$diagnosis <- relevel(factor(ss$diagnosis),"Control")
for (v in c("age","pmi","ph","prua")) {
  x <- suppressWarnings(as.numeric(as.character(ss[[v]]))); x[is.na(x)] <- median(x,na.rm=TRUE)
  ss[[paste0(v,"_z")]] <- as.numeric(scale(x))
}
cnt <- cnt[rowSums(cnt>=10) >= min(table(ss$diagnosis)),]
logf("GSE80655 unique donors: n=%d (%s) | %d genes", nrow(ss),
     paste(names(table(ss$diagnosis)), table(ss$diagnosis), sep="=", collapse=" "), nrow(cnt))

sym <- AnnotationDbi::mapIds(org.Hs.eg.db, sub("\\..*","",rownames(cnt)), "SYMBOL","ENSEMBL", multiVals="first")
ieg_id <- rownames(cnt)[!is.na(sym) & sym %in% IEG]
logf("IEGs present: %s", paste(sort(na.omit(sym[!is.na(sym) & sym %in% IEG])), collapse=", "))

## --- is pH itself different by diagnosis, and does it predict IEG level? ---
logf("\n--- pH and PMI by diagnosis ---")
for (v in c("ph","pmi")) {
  d <- data.frame(x=as.numeric(as.character(ss[[v]])), g=ss$diagnosis); d <- d[!is.na(d$x),]
  logf("  %-3s Ctrl %.2f | SCZ %.2f  Wilcoxon p = %.4f", v,
       median(d$x[d$g=="Control"]), median(d$x[d$g=="SCZ"]),
       wilcox.test(x~g,data=d)$p.value)
}
vsd <- assay(vst(DESeqDataSetFromMatrix(cnt, ss, ~1), blind=TRUE))
z <- t(scale(t(vsd[ieg_id,,drop=FALSE])))   # z-score each gene across samples
score <- colMeans(z, na.rm=TRUE)            # one value per sample
stopifnot(length(score)==nrow(ss))
logf("\n--- IEG module score ---")
logf("  vs pH   : r = %+.3f  p = %.4f", cor(score, ss$ph_z), cor.test(score, ss$ph_z)$p.value)
logf("  vs PMI  : r = %+.3f  p = %.4f", cor(score, ss$pmi_z), cor.test(score, ss$pmi_z)$p.value)
logf("  vs RIN-proxy PRUA: r = %+.3f", cor(score, ss$prua_z))
logf("  by diagnosis: Ctrl %.3f vs SCZ %.3f  Wilcoxon p = %.4f",
     median(score[ss$diagnosis=="Control"]), median(score[ss$diagnosis=="SCZ"]),
     wilcox.test(score ~ ss$diagnosis)$p.value)
logf("  partial (diagnosis | pH): p = %.4f",
     summary(lm(score ~ ph_z + pmi_z + diagnosis, data=ss))$coefficients["diagnosisSCZ",4])

## --- per-gene: with and without pH in the model ---
run <- function(form,lab){
  dds <- suppressMessages(DESeq(DESeqDataSetFromMatrix(cnt, ss, form), quiet=TRUE))
  r <- as.data.frame(results(dds, name=tail(resultsNames(dds),1)))
  r$sym <- sym[match(sub("\\..*","",rownames(r)), sub("\\..*","",rownames(cnt)))]
  s <- r[rownames(r) %in% ieg_id,]
  logf("\n  %s  (IEGs with nominal p<0.05: %d/%d)", lab, sum(s$pvalue<0.05,na.rm=TRUE), nrow(s))
  print(s[order(s$pvalue), c("sym","log2FoldChange","pvalue","padj")])
  invisible(s)
}
a <- run(~ age_z + pmi_z + prua_z + diagnosis,          "WITHOUT pH")
b <- run(~ age_z + pmi_z + ph_z + prua_z + diagnosis, "WITH pH")
pr <- read.csv(file.path(CS,"Deconvolution","cell_type_proportions_long.csv"), stringsAsFactors=FALSE)
pr <- subset(pr, signature=="MB" & cell_type %in% c("Neurons","Excitatory"))
pr <- aggregate(proportion ~ srr, pr, mean)
rt2 <- read.csv(RT, check.names=FALSE); mp2 <- setNames(rt2$BioSample, rt2$Run)
pr$bio <- ifelse(is.na(mp2[pr$srr]), pr$srr, mp2[pr$srr])
pr <- aggregate(proportion ~ bio, pr, mean)
ss$neuron_z <- as.numeric(scale(pr$proportion[match(rownames(ss), pr$bio)]))
logf("\nneuron proportion matched for %d/%d samples", sum(!is.na(ss$neuron_z)), nrow(ss))
if (!any(is.na(ss$neuron_z))) d2 <- run(~ age_z + pmi_z + ph_z + prua_z + neuron_z + diagnosis, "WITH pH + NEURONAL PROPORTION")
logf("\nmedian |log2FC| shrinkage after pH adjustment: %.1f%%",
     100*(1 - median(abs(b$log2FoldChange))/median(abs(a$log2FoldChange))))
