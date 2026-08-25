## GSE138082 with correct donor structure. CA1/CA3 only (DG metadata unreliable).
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
suppressMessages({library(DESeq2); library(limma); library(edgeR)})
logf <- function(...) cat(sprintf(...),"\n")
NP  <- "/rds/projects/e/elhamsak-dissertation-neuropsychiatric/RNA_seq"
CS  <- file.path(NP,"counts_scz"); SHD <- file.path(CS,"sheets")
OR  <- "/rds/projects/e/elhamsak-dissertation-autism-2/RNA_seq/DE_analysis"
OUT <- file.path(CS,"HIPP_donor"); dir.create(OUT, showWarnings=FALSE)

## ---- metadata straight from the series matrix ----
M <- file.path(SHD,"GSE138082_series_matrix.txt.gz")
ln <- readLines(gzfile(M)); rows <- list()
for (l in ln[grepl("^!Sample_", ln)]) {
  p <- gsub('"','', strsplit(l,"\t")[[1]]); k <- sub("^!","",p[1])
  if (k=="Sample_title") rows$title <- p[-1]
  else if (k=="Sample_characteristics_ch1") {
    tag <- trimws(sub(":.*","",p[2]))
    rows[[tag]] <- trimws(sub("^[^:]*:","",p[-1]))
  }
}
md <- as.data.frame(rows, stringsAsFactors=FALSE)
md$donor <- paste(md$gender, md$age, md$pmi, md$group, sep="_")

## sample sheet gives the SRR <-> title link via order in the matrix
ss <- read.delim(file.path(SHD,"GSE138082_samplesheet.tsv"), check.names=FALSE)
pr <- read.delim(file.path(SHD,"GSE138082_prua.tsv"))
stopifnot(nrow(ss)==nrow(md))
## join SRR -> GSM via the SRA run table, then GSM -> series-matrix row
RT <- "/rds/projects/e/elhamsak-dissertation-neuropsychiatric/RNA_seq/metadata/PRJNA574470_SraRunTable.csv"
rt <- read.csv(RT, check.names=FALSE)
gcol <- grep("GEO_Accession|Sample.Name|Sample Name", names(rt), value=TRUE)[1]
s2g  <- setNames(rt[[gcol]], rt$Run)
gsm  <- gsub('"','', strsplit(grep("^!Sample_geo_accession", ln, value=TRUE)[1], "\t")[[1]][-1])
md$gsm <- gsm
ss$gsm <- s2g[ss$srr]
logf("SRR->GSM mapped: %d/%d | matching series rows: %d",
     sum(!is.na(ss$gsm)), nrow(ss), sum(ss$gsm %in% md$gsm))
ss <- merge(ss, md[,c("gsm","title","hippocampal_subfield","drug","donor","diagnosis")],
            by="gsm", all.x=FALSE)
ss$group <- ifelse(grepl("_SZ_", ss$title), "SCZ", "Control")
rownames(ss) <- ss$srr
stopifnot(nrow(ss) > 0, length(unique(ss$group))==2)
print(table(ss$group, ss$hippocampal_subfield))
ss <- merge(ss, pr[,c("srr","prua")], by="srr")
rownames(ss) <- ss$srr
logf("loaded %d samples", nrow(ss))

## ---- counts ----
bt <- read.delim(file.path(OR,"gene_id2biotype.tsv"), header=FALSE, col.names=c("gid","bt"))
cm <- read.delim(file.path(NP,"PRJNA574470/counts/PRJNA574470_gene_counts_matrix.tsv"),
                 comment.char="#", check.names=FALSE)
g <- cm[[1]]; cm[[1]] <- NULL
cm <- cm[, !colnames(cm) %in% c("Chr","Start","End","Strand","Length"), drop=FALSE]
colnames(cm) <- regmatches(colnames(cm), regexpr("SRR[0-9]+", colnames(cm)))
cnt <- as.matrix(cm); rownames(cnt) <- g; storage.mode(cnt) <- "integer"
k <- intersect(colnames(cnt), rownames(ss)); cnt <- cnt[,k]; ss <- ss[k,]
cnt <- cnt[rownames(cnt) %in% bt$gid[bt$bt %in% c("protein_coding","lncRNA")],]

## ---- drop DG (metadata unreliable) ----
keep <- ss$hippocampal_subfield %in% c("CA1","CA3")
ss <- ss[keep,]; cnt <- cnt[,rownames(ss)]
logf("CA1/CA3 only: %d samples, %d donors (%s)", nrow(ss), length(unique(ss$donor)),
     paste(names(table(ss$group)), table(ss$group), sep="=", collapse=" "))
print(table(ss$group, ss$drug))

prep <- function(s,c){
  s$group <- relevel(factor(s$group),"Control")
  s$subfield <- factor(s$hippocampal_subfield); s$drug <- factor(s$drug)
  for (v in c("age","pmi","rin","prua")) s[[paste0(v,"_z")]] <- as.numeric(scale(as.numeric(s[[v]])))
  c <- c[rowSums(c>=10) >= min(table(s$group)),,drop=FALSE]
  list(s=s,c=c)
}
tally <- function(r,lab,lfc="log2FoldChange",p="padj"){
  s <- function(l) sum(!is.na(r[[p]]) & r[[p]]<0.05 & abs(r[[lfc]])>=l)
  logf("  %-42s padj05:%5d | >=0.3:%5d | >=1.0:%4d", lab, s(0), s(0.3), s(1))
  rownames(r)[!is.na(r[[p]]) & r[[p]]<0.05]
}

## ================= A. CA1 only, one sample per donor =================
logf("\n=========== A. CA1 ONLY (independent samples) ===========")
a <- prep(ss[ss$hippocampal_subfield=="CA1",], cnt[,ss$hippocampal_subfield=="CA1"])
logf("n=%d (%s) | %d genes", nrow(a$s),
     paste(names(table(a$s$group)), table(a$s$group), sep="=", collapse=" "), nrow(a$c))
fitD <- function(d,form,lab,tag){
  dds <- suppressMessages(DESeq(DESeqDataSetFromMatrix(d$c, d$s, form), quiet=TRUE))
  r <- as.data.frame(results(dds, name=tail(resultsNames(dds),1)))
  write.csv(r, file.path(OUT,paste0(tag,".csv")))
  tally(r,lab)
}
A1 <- fitD(a, ~ age_z+sex+pmi_z+rin_z+group,          "A1 age+sex+PMI+RIN",        "CA1_A1")
A2 <- fitD(a, ~ age_z+sex+pmi_z+rin_z+drug+group,     "A2 + drug covariate",       "CA1_A2")
off <- a$s$drug=="OFF"
ao <- prep(a$s[off,], a$c[,off])
logf("  drug-OFF subset: n=%d (%s)", nrow(ao$s),
     paste(names(table(ao$s$group)), table(ao$s$group), sep="=", collapse=" "))
A3 <- fitD(ao, ~ age_z+sex+pmi_z+rin_z+group,         "A3 drug-OFF only [SENSITIVITY]","CA1_A3")

## ================= B. CA1+CA3, donor-blocked =================
logf("\n=========== B. CA1+CA3, limma-voom + duplicateCorrelation ===========")
b <- prep(ss, cnt)
d <- calcNormFactors(DGEList(b$c, group=b$s$group))
des <- model.matrix(~ age_z+sex+pmi_z+rin_z+subfield+group, b$s)
v  <- voom(d, des)
cf <- duplicateCorrelation(v, des, block=b$s$donor)
logf("  within-donor correlation: %.3f", cf$consensus.correlation)
v  <- voom(d, des, block=b$s$donor, correlation=cf$consensus.correlation)
cf <- duplicateCorrelation(v, des, block=b$s$donor)
fit <- eBayes(lmFit(v, des, block=b$s$donor, correlation=cf$consensus.correlation))
tt  <- topTable(fit, coef="groupSCZ", number=Inf)
write.csv(tt, file.path(OUT,"CA1CA3_blocked.csv"))
B1 <- tally(tt, "B1 donor-blocked (+subfield)", "logFC", "adj.P.Val")

des2 <- model.matrix(~ age_z+sex+pmi_z+rin_z+subfield+drug+group, b$s)
v2 <- voom(d, des2, block=b$s$donor, correlation=cf$consensus.correlation)
cf2 <- duplicateCorrelation(v2, des2, block=b$s$donor)
fit2 <- eBayes(lmFit(v2, des2, block=b$s$donor, correlation=cf2$consensus.correlation))
tt2 <- topTable(fit2, coef="groupSCZ", number=Inf)
write.csv(tt2, file.path(OUT,"CA1CA3_blocked_drug.csv"))
B2 <- tally(tt2, "B2 donor-blocked + drug", "logFC", "adj.P.Val")

## ---- what happens if you IGNORE donor structure (the old analysis) ----
logf("\n--- for comparison: same samples, donor structure IGNORED ---")
fitN <- eBayes(lmFit(voom(d,des), des))
tally(topTable(fitN, coef="groupSCZ", number=Inf), "NAIVE (no blocking)", "logFC","adj.P.Val")

logf("\n=========== OVERLAP ===========")
logf("  A1 (CA1) n=%d | B1 (blocked) n=%d | intersect=%d",
     length(A1), length(B1), length(intersect(A1,B1)))
logf("  A1 vs A2(drug) intersect=%d | B1 vs B2(drug) intersect=%d",
     length(intersect(A1,A2)), length(intersect(B1,B2)))
logf("\nOutputs -> %s", OUT)

## ---- decomposition: is it DG, or is it the method? ----
logf("\n=========== DECOMPOSITION (DESeq2 throughout) ===========")
allss <- ss   # CA1+CA3 currently; rebuild full set
fitq <- function(s, c, lab){
  d <- prep(s,c)
  dds <- suppressMessages(DESeq(DESeqDataSetFromMatrix(d$c, d$s,
          ~ age_z+sex+pmi_z+rin_z+group), quiet=TRUE))
  r <- as.data.frame(results(dds, name=tail(resultsNames(dds),1)))
  tally(r, lab)
}
