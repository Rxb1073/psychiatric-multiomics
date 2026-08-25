#!/usr/bin/env Rscript
# Per-cohort composition-adjusted DE: does the IEG core survive astrocyte
# adjustment WITHIN each frontal cohort alone? Closes the LODO/composition gap.
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
suppressMessages(library(DESeq2))
logf <- function(...) cat(sprintf(...), "\n"); CASE <- "SCZ"
NEURO <- "/rds/projects/e/elhamsak-dissertation-neuropsychiatric/RNA_seq"
CS <- file.path(NEURO,"counts_scz"); SH <- file.path(CS,"sheets"); DEC <- file.path(CS,"Deconvolution")
OR <- "/rds/projects/e/elhamsak-dissertation-autism-2/RNA_seq/DE_analysis"
OUT <- file.path(CS,"PERCOHORT_composition"); dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
SIGNATURE <- "MB"; if (nzchar(Sys.getenv("DECONV_SIG"))) SIGNATURE <- Sys.getenv("DECONV_SIG")
logf("=== PER-COHORT COMPOSITION-ADJUSTED DE | signature: %s ===", SIGNATURE)
bt <- read.delim(file.path(OR,"gene_id2biotype.tsv"), header=FALSE, col.names=c("gid","bt")); btv <- setNames(bt$bt,bt$gid)
len <- read.delim(file.path(OR,"gene_lengths.tsv"), header=FALSE, col.names=c("gid","len")); lenv <- setNames(as.numeric(len$len),len$gid)
gm <- read.delim(file.path(OR,"gene_id2name.tsv"), header=FALSE, col.names=c("gid","nm")); nmv <- setNames(gm$nm,gm$gid)
IEG <- c("ATF3","NR4A1","NR4A2","NR4A3","EGR1","EGR2","FOS","FOSB","NPAS4","BTG2","ARC","DUSP1","JUNB")
ieg_ids <- names(nmv)[nmv %in% IEG]
read_counts <- function(p){ cm<-read.delim(p,comment.char="#",check.names=FALSE); g<-cm[[1]];cm[[1]]<-NULL
  cm<-cm[,!(colnames(cm)%in%c("Chr","Start","End","Strand","Length")),drop=FALSE]
  h<-regmatches(colnames(cm),regexpr("SRR[0-9]+",colnames(cm))); cm<-cm[,nchar(h)>0,drop=FALSE]
  colnames(cm)<-regmatches(colnames(cm),regexpr("SRR[0-9]+",colnames(cm))); rownames(cm)<-g
  m<-as.matrix(cm);m[is.na(m)]<-0;storage.mode(m)<-"integer";m }
props <- read.csv(file.path(DEC,"cell_type_proportions_long.csv"), stringsAsFactors=FALSE)
props <- subset(props, signature==SIGNATURE)
cohorts <- list(
  GSE80655_DLPFC=list(c=file.path(NEURO,"PRJNA319583/counts/PRJNA319583_gene_counts_matrix.tsv"),
    s=file.path(SH,"GSE80655_DLPFC_samplesheet.tsv"), base="~ age_z+sex+pmi_z+ph_z+diagnosis", num=c("age","pmi","ph")),
  GSE78936=list(c=file.path(NEURO,"PRJNA314463/counts/PRJNA314463_gene_counts.txt"),
    s=file.path(SH,"GSE78936_samplesheet.tsv"), base="~ brain_region+diagnosis", num=character(0)))
for (nm in names(cohorts)) {
  C <- cohorts[[nm]]; logf("\n############ %s ############", nm)
  ss <- read.delim(C$s, check.names=FALSE); rownames(ss) <- ss$srr
  cm <- read_counts(C$c); k <- intersect(colnames(cm),rownames(ss))
  cm <- cm[,k,drop=FALSE]; ss <- ss[k,,drop=FALSE]
  g <- rownames(cm)[!is.na(btv[rownames(cm)]) & btv[rownames(cm)]%in%c("protein_coding","lncRNA")]
  cm <- cm[g,,drop=FALSE]
  lnc <- g[btv[g]=="lncRNA"]; pc <- g[btv[g]=="protein_coding"]
  ss$q_lnc <- colSums(cm[lnc,,drop=FALSE])/colSums(cm)
  Lp <- lenv[pc]; ok <- !is.na(Lp)&Lp>0; W <- cm[pc,,drop=FALSE][ok,,drop=FALSE]; Lo <- Lp[ok]
  ss$q_len <- apply(W,2,function(x){s<-sum(x); if(is.na(s)||s==0) NA_real_ else sum(x*Lo)/s})
  pw <- subset(props, dataset==nm)
  if (nrow(pw)==0) pw <- subset(props, dataset==sub("_DLPFC","",nm))
  w <- reshape(pw[,c("srr","cell_type","proportion")], idvar="srr", timevar="cell_type", direction="wide")
  names(w) <- sub("^proportion\\.","",names(w)); rownames(w) <- w$srr
  ss$astro_z <- as.numeric(scale(w[rownames(ss),"Astrocytes"]))
  ncol_ <- if ("Neurons" %in% names(w)) "Neurons" else "Excitatory"
  ss$neuron_z <- as.numeric(scale(w[rownames(ss),ncol_]))
  ss$q_len_z <- as.numeric(scale(ss$q_len)); ss$q_lnc_z <- as.numeric(scale(ss$q_lnc))
  for (v in C$num) ss[[paste0(v,"_z")]] <- as.numeric(scale(suppressWarnings(as.numeric(ss[[v]]))))
  ss$diagnosis <- relevel(factor(ss$diagnosis),"Control")
  if ("sex" %in% names(ss)) ss$sex <- factor(ss$sex)
  if ("brain_region" %in% names(ss)) ss$brain_region <- droplevels(factor(ss$brain_region))
  vars <- all.vars(as.formula(C$base))
  vars <- c(vars,"astro_z","neuron_z","q_len_z","q_lnc_z")
  vars <- intersect(vars, names(ss))
  keep <- complete.cases(ss[,vars,drop=FALSE])
  if (any(!keep)) logf("  dropping %d samples with NA", sum(!keep))
  ss <- ss[keep,]; cm <- cm[,rownames(ss),drop=FALSE]
  sg <- min(table(ss$diagnosis)); cm <- cm[rowSums(cm>=10)>=sg,,drop=FALSE]
  logf("  n=%d (%s) | %d genes | cor(astro_z,q_len_z)=%+.3f", ncol(cm),
       paste(names(table(ss$diagnosis)),table(ss$diagnosis),sep="=",collapse="/"), nrow(cm),
       cor(ss$astro_z, ss$q_len_z, use="complete.obs"))
  logf("  astro_z by dx: Control %.3f | SCZ %.3f (Wilcoxon p=%.4f)",
       median(ss$astro_z[ss$diagnosis=="Control"]), median(ss$astro_z[ss$diagnosis=="SCZ"]),
       suppressWarnings(wilcox.test(astro_z~diagnosis,data=ss)$p.value))
  fit <- function(f,lab){
    dds <- suppressMessages(DESeq(DESeqDataSetFromMatrix(cm, ss, as.formula(f)), quiet=TRUE))
    r <- results(dds, contrast=c("diagnosis",CASE,"Control"), alpha=0.05)
    s <- subset(as.data.frame(r), !is.na(padj)&padj<0.05)
    logf("  %-46s : %4d genes (up=%d down=%d)", lab, nrow(s), sum(s$log2FoldChange>0), sum(s$log2FoldChange<0))
    list(res=r, sig=rownames(s)) }
  b <- C$base
  c0 <- fit(b, "C0 baseline (own covariates)")
  c1 <- fit(paste(b,"+ astro_z"), "C1 + astrocyte [MEDIATION]")
  c2 <- fit(paste(b,"+ q_len_z + q_lnc_z"), "C2 + quality")
  c3 <- fit(paste(b,"+ astro_z + q_len_z + q_lnc_z"), "C3 + both")
  logf("  fate of C0's %d genes: C1 %d | C2 %d | C3 %d", length(c0$sig),
       length(intersect(c0$sig,c1$sig)), length(intersect(c0$sig,c2$sig)), length(intersect(c0$sig,c3$sig)))
  logf("\n  --- IEG fate in %s (padj per model) ---", nm)
  hdr <- sprintf("  %-8s %8s %8s %8s %8s %9s", "gene","C0","C1","C2","C3","LFC(C3)")
  logf("%s", hdr)
  for (id in ieg_ids) {
    if (!id %in% rownames(c0$res)) next
    p <- function(m) { v <- as.data.frame(m$res)[id,"padj"]; if (is.na(v)) "  NA" else sprintf("%.4f", v) }
    l3 <- as.data.frame(c3$res)[id,"log2FoldChange"]
    logf("  %-8s %8s %8s %8s %8s %9.3f", nmv[id], p(c0), p(c1), p(c2), p(c3), l3) }
  saveRDS(list(c0=c0,c1=c1,c2=c2,c3=c3,meta=ss), file.path(OUT, sprintf("%s_%s.rds", nm, SIGNATURE)))
}
logf("\n=== READ: if IEG padj stays <0.05 under C1/C3, the IEG signal is cell-intrinsic")
logf("    WITHIN that cohort alone - i.e. survives composition adjustment AND is not pooled-only.")
sink(file.path(OUT,"sessionInfo.txt")); print(sessionInfo()); sink()
logf("\noutputs -> %s", OUT)
