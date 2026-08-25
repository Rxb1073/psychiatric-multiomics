#!/usr/bin/env Rscript
## Apply each ORIGINAL STUDY'S OWN methodology to donor-correct samples.
## Output -> counts_scz/ORIGINAL_METHODS/  (nothing existing is overwritten)
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
suppressMessages({
  library(DESeq2); library(edgeR); library(limma); library(sva)
  library(AnnotationDbi); library(org.Hs.eg.db)
})
logf <- function(...) cat(sprintf(...), "\n")
NP  <- "/rds/projects/e/elhamsak-dissertation-neuropsychiatric/RNA_seq"
CS  <- file.path(NP, "counts_scz"); SHD <- file.path(CS, "sheets")
OR  <- "/rds/projects/e/elhamsak-dissertation-autism-2/RNA_seq/DE_analysis"
RT8 <- "/rds/projects/e/elhamsak-dissertation-autism-2/MDD/transcriptomics/GSE80655/metadata/PRJNA319583_SraRunTable.csv"
OUT <- file.path(CS, "ORIGINAL_METHODS"); dir.create(OUT, showWarnings = FALSE)
bt   <- read.delim(file.path(OR, "gene_id2biotype.tsv"), header = FALSE, col.names = c("gid","bt"))
KEEP <- bt$gid[bt$bt %in% c("protein_coding","lncRNA")]
read_counts <- function(p) {
  cm <- read.delim(p, comment.char = "#", check.names = FALSE)
  g  <- cm[[1]]; cm[[1]] <- NULL
  cm <- cm[, !colnames(cm) %in% c("Chr","Start","End","Strand","Length"), drop = FALSE]
  colnames(cm) <- regmatches(colnames(cm), regexpr("SRR[0-9]+", colnames(cm)))
  m <- as.matrix(cm); rownames(m) <- g; storage.mode(m) <- "integer"; m
}
sym_of <- function(ids) {
  s <- suppressMessages(AnnotationDbi::mapIds(org.Hs.eg.db, sub("\\..*","",ids),
                        "SYMBOL","ENSEMBL", multiVals="first"))
  ifelse(is.na(s), ids, s)
}
report <- function(res, lab, tag, pcol="padj", fcol="log2FoldChange") {
  res$symbol <- sym_of(rownames(res))
  write.csv(res, file.path(OUT, paste0(tag,"_full.csv")))
  sig <- res[!is.na(res[[pcol]]) & res[[pcol]] < 0.05, ]
  sig <- sig[order(sig[[pcol]]), ]
  write.csv(sig, file.path(OUT, paste0(tag,"_sig.csv")))
  logf("  %-46s FDR<0.05: %4d   (up %d / down %d)", lab, nrow(sig),
       sum(sig[[fcol]] > 0, na.rm=TRUE), sum(sig[[fcol]] < 0, na.rm=TRUE))
  if (nrow(sig)) print(utils::head(sig[, c("symbol", fcol, pcol)], 15))
  rownames(sig)
}

logf("\n============ CA1 (GSE138082) : log-CPM + SVA + lm ============")
ln <- readLines(gzfile(file.path(SHD,"GSE138082_series_matrix.txt.gz"))); rows <- list()
for (l in ln[grepl("^!Sample_", ln)]) {
  p <- gsub('"','', strsplit(l,"\t")[[1]]); k <- sub("^!","",p[1])
  if (k=="Sample_title") rows$title <- p[-1]
  else if (k=="Sample_characteristics_ch1") {
    tag <- trimws(sub(":.*","",p[2])); rows[[tag]] <- trimws(sub("^[^:]*:","",p[-1])) }
}
md <- as.data.frame(rows, stringsAsFactors=FALSE)
md$gsm <- gsub('"','', strsplit(grep("^!Sample_geo_accession", ln, value=TRUE)[1],"\t")[[1]][-1])
rt   <- read.csv(file.path(NP,"metadata/PRJNA574470_SraRunTable.csv"), check.names=FALSE)
gcol <- grep("GEO_Accession|Sample.Name|Sample Name", names(rt), value=TRUE)[1]
hs   <- read.delim(file.path(SHD,"GSE138082_samplesheet.tsv"), check.names=FALSE)
hs$gsm <- setNames(rt[[gcol]], rt$Run)[hs$srr]
hs <- merge(hs, md[,c("gsm","title","hippocampal_subfield","drug")], by="gsm")
hs$group <- ifelse(grepl("_SZ_", hs$title), "SCZ", "Control")
rownames(hs) <- hs$srr; hs <- hs[hs$hippocampal_subfield=="CA1",]
hc <- read_counts(file.path(NP,"PRJNA574470/counts/PRJNA574470_gene_counts_matrix.tsv"))
k <- intersect(colnames(hc), rownames(hs)); hc <- hc[,k]; hs <- hs[k,]
hc <- hc[rownames(hc) %in% KEEP,]
cpm0 <- cpm(DGEList(hc))
keep <- apply(sapply(unique(hs$group), function(g)
          rowSums(cpm0[, hs$group==g, drop=FALSE] > 0) == sum(hs$group==g)), 1, any)
hc <- hc[keep,]
logf("n = %d (%s) | %d genes after source CPM filter", nrow(hs),
     paste(names(table(hs$group)), table(hs$group), sep="=", collapse=" "), nrow(hc))
hs$group <- relevel(factor(hs$group), "Control")
for (v in c("age","pmi","rin")) hs[[v]] <- as.numeric(as.character(hs[[v]]))
E <- log2(cpm(DGEList(hc)) + 1)
mod  <- model.matrix(~ group + age + sex + pmi + rin, hs)
mod0 <- model.matrix(~ age + sex + pmi + rin, hs)
nsv  <- tryCatch(num.sv(E, mod, method="be", B=100), error=function(e) 0)
logf("SVA: %d surrogate variables estimated", nsv)
des <- mod
if (nsv > 0) {
  sv <- sva(E, mod, mod0, n.sv=nsv)$sv
  des <- cbind(mod, sv); colnames(des)[(ncol(mod)+1):ncol(des)] <- paste0("SV",1:nsv)
}
ca1_sva <- report(topTable(eBayes(lmFit(E, des)), coef="groupSCZ", number=Inf, sort.by="none"),
  sprintf("CA1 lm + %d SVs [GSE138082 method]", nsv), "CA1_GSE138082method", "adj.P.Val","logFC")
report(topTable(eBayes(lmFit(E, mod)), coef="groupSCZ", number=Inf, sort.by="none"),
  "CA1 lm, no SVs [reference]", "CA1_noSVA", "adj.P.Val","logFC")

logf("\n============ DLPFC (GSE80655) : DESeq2 LRT, 4 covariates ============")
d8 <- read.delim(file.path(SHD,"GSE80655_DLPFC_samplesheet.tsv"), check.names=FALSE)
d8 <- merge(d8, read.delim(file.path(SHD,"GSE80655_DLPFC_prua.tsv"))[,c("srr","prua")], by="srr")
rownames(d8) <- d8$srr
c8 <- read_counts(file.path(NP,"PRJNA319583/counts/PRJNA319583_gene_counts_matrix.tsv"))
k <- intersect(colnames(c8), rownames(d8)); c8 <- c8[,k]; d8 <- d8[k,]
rt8 <- read.csv(RT8, check.names=FALSE); bio <- setNames(rt8$BioSample, rt8$Run)[colnames(c8)]
idx <- split(seq_len(ncol(c8)), bio)
c8s <- vapply(idx, function(i) as.integer(rowSums(c8[,i,drop=FALSE])), integer(nrow(c8)))
rownames(c8s) <- rownames(c8); ke <- !duplicated(bio)
d8 <- d8[ke,]; rownames(d8) <- bio[ke]; c8 <- c8s[, rownames(d8)]
c8 <- c8[rownames(c8) %in% KEEP,]
d8$diagnosis <- relevel(factor(d8$diagnosis), "Control")
for (v in c("age","pmi","ph","prua")) {
  x <- suppressWarnings(as.numeric(as.character(d8[[v]])))
  x[is.na(x)] <- median(x, na.rm=TRUE)
  d8[[v]] <- as.numeric(scale(x))          # centre and scale: avoids non-convergence
}
logf("n = %d (%s) | %d genes", nrow(d8),
     paste(names(table(d8$diagnosis)), table(d8$diagnosis), sep="=", collapse=" "), nrow(c8))
dds <- suppressMessages(DESeq(DESeqDataSetFromMatrix(c8, d8, ~ age+pmi+ph+prua+diagnosis),
        test="LRT", reduced = ~ age+pmi+ph+prua, quiet=TRUE))
conv <- mcols(dds)$fullBetaConv
if (is.null(conv)) conv <- rep(TRUE, nrow(dds))
conv[is.na(conv)] <- FALSE
logf("non-convergent genes removed (as in source): %d", sum(!conv))
dlpfc_lrt <- report(as.data.frame(results(dds)[which(conv),]),
  "DLPFC DESeq2 LRT [GSE80655 method]", "DLPFC_GSE80655method")

logf("\n============ BA11 (GSE78936) : union of three methods ============")
d7 <- read.delim(file.path(SHD,"GSE78936_samplesheet.tsv"), check.names=FALSE)
d7 <- merge(d7, read.delim(file.path(SHD,"GSE78936_prua.tsv"))[,c("srr","prua")], by="srr")
rownames(d7) <- d7$srr; d7 <- d7[d7$brain_region=="Brodmann's Area 11",]
c7 <- read_counts(file.path(NP,"PRJNA314463/counts/PRJNA314463_gene_counts.txt"))
k <- intersect(colnames(c7), rownames(d7)); c7 <- c7[,k]; d7 <- d7[k,]
c7 <- c7[rownames(c7) %in% KEEP,]
c7 <- c7[rowSums(c7 >= 2) >= 0.5*ncol(c7), ]
d7$diagnosis <- relevel(factor(d7$diagnosis), "Control")
logf("n = %d (%s) | %d genes after source filter", nrow(d7),
     paste(names(table(d7$diagnosis)), table(d7$diagnosis), sep="=", collapse=" "), nrow(c7))
des7 <- model.matrix(~ diagnosis, d7)
s1 <- report(as.data.frame(results(suppressMessages(
        DESeq(DESeqDataSetFromMatrix(c7, d7, ~ diagnosis), quiet=TRUE)))),
        "  (a) DESeq2", "BA11_DESeq2")
y  <- estimateDisp(calcNormFactors(DGEList(c7, group=d7$diagnosis)), des7)
r2 <- topTags(glmQLFTest(glmQLFit(y, des7), coef=2), n=Inf)$table
colnames(r2)[colnames(r2)=="FDR"] <- "padj"; colnames(r2)[colnames(r2)=="logFC"] <- "log2FoldChange"
s2 <- report(r2, "  (b) edgeR", "BA11_edgeR")
s3 <- report(topTable(eBayes(lmFit(voom(y, des7), des7)), coef=2, number=Inf, sort.by="none"),
        "  (c) voom-limma", "BA11_voomlimma", "adj.P.Val","logFC")
uni <- unique(c(s1,s2,s3))
logf("\n  UNION (source rule, sig in >=1 method): %d genes", length(uni))
logf("    DESeq2 %d | edgeR %d | voom %d | all three: %d",
     length(s1), length(s2), length(s3), length(Reduce(intersect, list(s1,s2,s3))))
if (length(uni)) {
  u <- data.frame(gene=uni, symbol=sym_of(uni),
                  DESeq2=uni %in% s1, edgeR=uni %in% s2, voom=uni %in% s3)
  write.csv(u, file.path(OUT,"BA11_GSE78936method_union.csv"), row.names=FALSE)
  print(utils::head(u, 20))
}

logf("\n============ SUMMARY (source methods, donor-correct, FDR<0.05) ============")
logf("  CA1   GSE138082 lm + SVA      : %4d genes  [26 donors]", length(ca1_sva))
logf("  DLPFC GSE80655 DESeq2 LRT     : %4d genes  [48 donors]", length(dlpfc_lrt))
logf("  BA11  GSE78936 union of three : %4d genes  [28 donors]", length(uni))
ov <- Reduce(intersect, list(sub("\\..*","",ca1_sva), sub("\\..*","",dlpfc_lrt), sub("\\..*","",uni)))
logf("  shared by all three           : %4d genes", length(ov))
if (length(ov)) print(data.frame(gene=ov, symbol=sym_of(ov)))
logf("\nNo fold-change filter (none of the source studies used one).")
logf("Outputs -> %s", OUT)
