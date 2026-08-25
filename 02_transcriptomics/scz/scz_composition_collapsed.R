#!/usr/bin/env Rscript
# SCZ composition-adjusted DE. Adapted from ASD pooled_DE_composition_adjusted.R
# MEDIATION (astrocyte) vs CONFOUND (quality) — see ASD script's interpretive note.
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
suppressMessages(library(DESeq2))
logf <- function(...) cat(sprintf(...), "\n"); CASE <- "SCZ"
NEURO <- "/rds/projects/e/elhamsak-dissertation-neuropsychiatric/RNA_seq"
CS <- file.path(NEURO,"counts_scz"); SH <- file.path(CS,"sheets")
DEC <- file.path(CS,"Deconvolution")
OR <- "/rds/projects/e/elhamsak-dissertation-autism-2/RNA_seq/DE_analysis"
OUT <- file.path(CS,"POOLED_composition_COLLAPSED"); dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
SIGNATURE <- "MB"; if (nzchar(Sys.getenv("DECONV_SIG"))) SIGNATURE <- Sys.getenv("DECONV_SIG")
logf("using deconvolution signature: %s", SIGNATURE)
bt <- read.delim(file.path(OR,"gene_id2biotype.tsv"), header=FALSE, col.names=c("gid","bt")); btv <- setNames(bt$bt,bt$gid)
len <- read.delim(file.path(OR,"gene_lengths.tsv"), header=FALSE, col.names=c("gid","len")); lenv <- setNames(as.numeric(len$len),len$gid)
gm <- read.delim(file.path(OR,"gene_id2name.tsv"), header=FALSE, col.names=c("gid","nm")); nmv <- setNames(gm$nm,gm$gid)
read_counts <- function(p){ cm<-read.delim(p,comment.char="#",check.names=FALSE); g<-cm[[1]];cm[[1]]<-NULL
  cm<-cm[,!(colnames(cm)%in%c("Chr","Start","End","Strand","Length")),drop=FALSE]
  h<-regmatches(colnames(cm),regexpr("SRR[0-9]+",colnames(cm))); cm<-cm[,nchar(h)>0,drop=FALSE]
  colnames(cm)<-regmatches(colnames(cm),regexpr("SRR[0-9]+",colnames(cm))); rownames(cm)<-g
  m<-as.matrix(cm);m[is.na(m)]<-0;storage.mode(m)<-"integer";m }
ds <- list(
  list(c=file.path(NEURO,"PRJNA319583/counts/PRJNA319583_gene_counts_matrix.tsv"), s=file.path(SH,"GSE80655_DLPFC_samplesheet.tsv"), n="GSE80655"),
  list(c=file.path(NEURO,"PRJNA314463/counts/PRJNA314463_gene_counts.txt"), s=file.path(SH,"GSE78936_samplesheet.tsv"), n="GSE78936"))
M<-list(); T<-list(); MEMB<-list()
for (d in ds) {
  ss <- read.delim(d$s, check.names=FALSE); rownames(ss) <- ss$srr
  cm <- read_counts(d$c); k <- intersect(colnames(cm), rownames(ss))
  cm <- cm[,k,drop=FALSE]; ss <- ss[k,,drop=FALSE]
  RT <- "/rds/projects/e/elhamsak-dissertation-autism-2/MDD/transcriptomics/GSE80655/metadata/PRJNA319583_SraRunTable.csv"
  if (file.exists(RT)) {
    rt  <- read.csv(RT, check.names=FALSE); mp <- setNames(rt$BioSample, rt$Run)
    bio <- mp[colnames(cm)]
    if (!any(is.na(bio)) && length(unique(bio)) < ncol(cm)) {
      logf("  %s run-splits: %d runs -> %d samples (summing)", d$n, ncol(cm), length(unique(bio)))
      idx <- split(seq_len(ncol(cm)), bio)
      cm2 <- vapply(idx, function(i) as.integer(rowSums(cm[,i,drop=FALSE])), integer(nrow(cm)))
      rownames(cm2) <- rownames(cm)
      for (nm in names(idx)) MEMB[[nm]] <- colnames(cm)[idx[[nm]]]
      keep <- !duplicated(bio); ss <- ss[keep,,drop=FALSE]; ss$srr <- bio[keep]
      rownames(ss) <- ss$srr; cm <- cm2[, ss$srr, drop=FALSE]
      storage.mode(cm) <- "integer"
    }
  }
  M[[d$n]] <- cm
  T[[d$n]] <- data.frame(srr=ss$srr, diagnosis=ss$diagnosis, study=d$n,
    brain_region=ss$brain_region, stringsAsFactors=FALSE) }
g <- Reduce(intersect, lapply(M,rownames)); g <- g[!is.na(btv[g]) & btv[g]%in%c("protein_coding","lncRNA")]
cnt <- do.call(cbind, lapply(M, function(m) m[g,,drop=FALSE]))
mt <- do.call(rbind,T); rownames(mt) <- mt$srr; mt <- mt[colnames(cnt),]
logf("pooled frontal: %d samples x %d genes", ncol(cnt), nrow(cnt))
lnc <- g[btv[g]=="lncRNA"]; pc <- g[btv[g]=="protein_coding"]
mt$q_lnc <- colSums(cnt[lnc,,drop=FALSE])/colSums(cnt)
Lp <- lenv[pc]; ok <- !is.na(Lp)&Lp>0; W <- cnt[pc,,drop=FALSE][ok,,drop=FALSE]; Lo <- Lp[ok]
mt$q_len <- apply(W,2,function(x){s<-sum(x); if(is.na(s)||s==0) NA_real_ else sum(x*Lo)/s})
props <- read.csv(file.path(DEC,"cell_type_proportions_long.csv"), stringsAsFactors=FALSE)
props <- subset(props, signature==SIGNATURE)
wide <- reshape(props[,c("srr","cell_type","proportion")], idvar="srr", timevar="cell_type", direction="wide")
names(wide) <- sub("^proportion\\.","",names(wide)); rownames(wide) <- wide$srr
if (length(MEMB)) {
  num <- setdiff(names(wide),"srr")
  add <- do.call(rbind, lapply(names(MEMB), function(nm){
    r <- wide[intersect(MEMB[[nm]], rownames(wide)), num, drop=FALSE]
    if (!nrow(r)) return(NULL)
    d <- as.data.frame(t(colMeans(r, na.rm=TRUE))); d$srr <- nm; d[, c("srr",num)] }))
  if (!is.null(add)) { rownames(add) <- add$srr; wide <- rbind(wide[,c("srr",num)], add) }
  logf("  averaged deconvolution proportions for %d collapsed samples", length(MEMB))
}
if (!"Astrocytes" %in% names(wide)) stop("Astrocytes absent from deconv output")
mt$astro <- wide[rownames(mt),"Astrocytes"]
neuron_col <- if ("Neurons" %in% names(wide)) "Neurons" else "Excitatory"
mt$neuron <- wide[rownames(mt),neuron_col]; logf("neuronal column used: %s", neuron_col)
zwithin <- function(x,grp) ave(x,grp,FUN=function(v) as.numeric(scale(v)))
mt$astro_z <- zwithin(mt$astro,mt$study); mt$neuron_z <- zwithin(mt$neuron,mt$study)
mt$q_len_z <- zwithin(mt$q_len,mt$study); mt$q_lnc_z <- zwithin(mt$q_lnc,mt$study)
mt$diagnosis <- relevel(factor(mt$diagnosis),"Control"); mt$study <- factor(mt$study)
logf("\ncor(astro_z,q_len_z)=%+.3f | cor(astro_z,q_lnc_z)=%+.3f",
  cor(mt$astro_z,mt$q_len_z,use="complete.obs"), cor(mt$astro_z,mt$q_lnc_z,use="complete.obs"))
logf("astrocyte z by diagnosis: Control %.3f | SCZ %.3f (Wilcoxon p=%.4f)",
  median(mt$astro_z[mt$diagnosis=="Control"],na.rm=TRUE), median(mt$astro_z[mt$diagnosis=="SCZ"],na.rm=TRUE),
  suppressWarnings(wilcox.test(astro_z~diagnosis,data=mt)$p.value))
keep <- complete.cases(mt[,c("diagnosis","astro_z","neuron_z","q_len_z","q_lnc_z")])
cnt <- cnt[,keep,drop=FALSE]; mt <- mt[keep,]
sg <- min(table(mt$diagnosis)); cnt <- cnt[rowSums(cnt>=10)>=sg,,drop=FALSE]
logf("after filtering: %d samples x %d genes\n", ncol(cnt), nrow(cnt))
fit <- function(form,lab,meta=mt,counts=cnt){
  dds <- suppressMessages(DESeq(DESeqDataSetFromMatrix(counts,meta,form), quiet=TRUE))
  r <- results(dds, contrast=c("diagnosis",CASE,"Control"), alpha=0.05)
  s <- subset(as.data.frame(r), !is.na(padj)&padj<0.05)
  logf("%-56s : %4d genes (up=%d, down=%d)", lab, nrow(s), sum(s$log2FoldChange>0), sum(s$log2FoldChange<0))
  list(res=r, sig=rownames(s)) }
logf("=========== POOLED FRONTAL MODELS (n=%d) ===========", ncol(cnt))
m0 <- fit(~ study + diagnosis,                                  "M0 baseline (~study+diagnosis)")
m1 <- fit(~ study + astro_z + diagnosis,                        "M1 + astrocyte  [MEDIATION]")
m2 <- fit(~ study + q_len_z + q_lnc_z + diagnosis,              "M2 + quality proxies [CONFOUND]")
m3 <- fit(~ study + astro_z + q_len_z + q_lnc_z + diagnosis,    "M3 + both")
m4 <- fit(~ study + astro_z + neuron_z + diagnosis,             "M4 + astrocyte + neuronal [COMPOSITION]")
base <- m0$sig
logf("\n=========== FATE OF THE %d BASELINE GENES ===========", length(base))
for (nm in c("M1 astro","M2 quality","M3 both","M4 composition")) {
  s <- switch(nm,"M1 astro"=m1$sig,"M2 quality"=m2$sig,"M3 both"=m3$sig,"M4 composition"=m4$sig)
  logf("  retained under %-16s : %4d / %d (%.0f%%)", nm, length(intersect(base,s)), length(base),
       100*length(intersect(base,s))/max(1,length(base))) }
surv <- Reduce(intersect, list(base,m1$sig,m2$sig,m3$sig,m4$sig))
logf("\n  survive EVERY model: %d", length(surv))
if (length(surv)) { r <- as.data.frame(m3$res)[surv,]
  r$gene <- ifelse(is.na(nmv[surv]), surv, nmv[surv]); r$biotype <- btv[surv]
  print(head(r[order(r$padj), c("gene","biotype","log2FoldChange","padj")], 30)) }
write.table(data.frame(gene_id=base, gene=nmv[base], biotype=btv[base],
  in_M1=base%in%m1$sig, in_M2=base%in%m2$sig, in_M3=base%in%m3$sig, in_M4=base%in%m4$sig),
  file.path(OUT,sprintf("baseline_fate_%s.tsv",SIGNATURE)), sep="\t", quote=FALSE, row.names=FALSE)
saveRDS(list(m0=m0,m1=m1,m2=m2,m3=m3,m4=m4,meta=mt), file.path(OUT,sprintf("models_%s.rds",SIGNATURE)))
# ---- HIPPOCAMPUS arm: real RIN available (analogue of ASD GSE102741) ----
logf("\n\n=========== GSE138082 HIPPOCAMPUS - REAL RIN ===========")
hs <- read.delim(file.path(SH,"GSE138082_samplesheet.tsv"), check.names=FALSE); rownames(hs) <- hs$srr
hc <- read_counts(file.path(NEURO,"PRJNA574470/counts/PRJNA574470_gene_counts_matrix.tsv"))
k <- intersect(colnames(hc), rownames(hs)); hc <- hc[,k,drop=FALSE]; hs <- hs[k,,drop=FALSE]
hg <- rownames(hc)[!is.na(btv[rownames(hc)]) & btv[rownames(hc)]%in%c("protein_coding","lncRNA")]
hc <- hc[hg,,drop=FALSE]
hp <- subset(props, dataset=="GSE138082")
hw <- reshape(hp[,c("srr","cell_type","proportion")], idvar="srr", timevar="cell_type", direction="wide")
names(hw) <- sub("^proportion\\.","",names(hw)); rownames(hw) <- hw$srr
hs$astro_z <- as.numeric(scale(hw[rownames(hs),"Astrocytes"]))
hs$diagnosis <- relevel(factor(hs$diagnosis),"Control"); hs$sex <- factor(hs$sex)
hs$age_z <- as.numeric(scale(suppressWarnings(as.numeric(hs$age))))
hs$pmi_z <- as.numeric(scale(suppressWarnings(as.numeric(hs$pmi))))
hs$rin_z <- as.numeric(scale(suppressWarnings(as.numeric(hs$rin))))
ok2 <- complete.cases(hs[,c("age_z","sex","pmi_z","rin_z","astro_z","diagnosis")])
hs <- hs[ok2,]; hc <- hc[,rownames(hs),drop=FALSE]
sg2 <- min(table(hs$diagnosis)); hc <- hc[rowSums(hc>=10)>=sg2,,drop=FALSE]
logf("n=%d samples, %d genes | cor(astro_z,rin)=%+.3f", ncol(hc), nrow(hc),
     cor(hs$astro_z, hs$rin_z, use="complete.obs"))
r0 <- fit(~ age_z+sex+pmi_z+diagnosis,                  "R0 age+sex+pmi",             hs, hc)
r1 <- fit(~ age_z+sex+pmi_z+rin_z+diagnosis,            "R1 + RIN (measured)",        hs, hc)
r2 <- fit(~ age_z+sex+pmi_z+rin_z+astro_z+diagnosis,    "R2 + RIN + astrocyte",       hs, hc)
logf("\n  hippocampus R0 -> R1 retained: %d/%d | R0 -> R2 retained: %d/%d",
     length(intersect(r0$sig,r1$sig)), length(r0$sig), length(intersect(r0$sig,r2$sig)), length(r0$sig))
logf("\n=========== HOW TO READ THIS ===========")
logf(" M1 (astrocyte) kills genes -> difference MEDIATED by composition (astrogliosis), not artefact.")
logf(" M2 (quality) kills them but M1 doesn't -> degradation CONFOUND (artefact).")
logf(" Both kill -> ambiguous (q_len/q_lnc correlate with composition; see printed cor above).")
logf(" Hippocampus R1 uses MEASURED RIN -> cleanest quality test; R2 adds composition.")
sink(file.path(OUT,"sessionInfo.txt")); print(sessionInfo()); sink()
logf("\noutputs -> %s", OUT)

## ================= PER-COHORT MEDIATION LADDERS =================
run_ladder <- function(sel, lab) {
  mm <- mt[sel, , drop=FALSE]; cc <- cnt[, rownames(mm), drop=FALSE]
  mm$diagnosis <- droplevels(factor(mm$diagnosis))
  if (length(unique(mm$diagnosis)) < 2 || min(table(mm$diagnosis)) < 5) {
    logf("  %s: too few samples, skipped", lab); return(invisible(NULL)) }
  for (v in c("astro","neuron","q_len","q_lnc"))
    mm[[paste0(v,"_z")]] <- as.numeric(scale(mm[[v]]))
  cc <- cc[rowSums(cc >= 10) >= min(table(mm$diagnosis)), , drop=FALSE]
  logf("\n--- %s: n=%d (%s) | %d genes ---", lab, nrow(mm),
       paste(names(table(mm$diagnosis)), table(mm$diagnosis), sep="=", collapse=" "), nrow(cc))
  fit(~ diagnosis,                                  paste(lab,"P0 baseline"),      mm, cc)
  fit(~ astro_z + diagnosis,                        paste(lab,"P1 +astro [MED]"),  mm, cc)
  fit(~ q_len_z + q_lnc_z + diagnosis,              paste(lab,"P2 +quality"),      mm, cc)
  fit(~ astro_z + neuron_z + diagnosis,             paste(lab,"P4 +composition"),  mm, cc)
}
logf("\n=========== PER-COHORT LADDERS (collapsed) ===========")
run_ladder(mt$study=="GSE80655", "GSE80655")
run_ladder(mt$study=="GSE78936" & mt$brain_region=="Brodmann's Area 11", "GSE78936-BA11")
